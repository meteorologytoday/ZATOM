using Formatting
using ArgMacros

parsed = @dictarguments begin
    @argumentrequired String  label       "--label"
    @argumentrequired String  p_name      "--p-name"
    @argumentrequired Float64 p_val       "--p-val"
    @argumentrequired Float64 p_fac       "--p-fac"
    @argumentrequired Int64   ds_exp_lb   "--ds-exp-lb"
    @argumentrequired String  scan_dir    "--scan-dir"
    @argumentrequired String  shape       "--shape"
    @argumentrequired Int64   ncpu        "--ncpu"
    @argumentrequired Int64   nthread     "--nthread"
end

scan_dir = parsed[:scan_dir]

if ! ( scan_dir in ["pos", "neg"] )
    throw(ErrorException("Unknown --scan-dir : $scan_dir"))
end

casename = format("{:s}{:s}_{:s}_P{:.2f}_lb{:d}_{:s}", parsed[:label], parsed[:shape], parsed[:p_name], parsed[:p_val] * parsed[:p_fac], parsed[:ds_exp_lb], scan_dir)
open("submit_$(casename).sh", "w") do io

    write(io,
"""
#!/bin/bash

#SBATCH -p cw3e-compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH -t 24:00:00
#SBATCH -J $casename
#SBATCH -A csg102
#SBATCH -o R-$(casename).%j.%N.out
#SBATCH -e R-$(casename).%j.%N.err
#SBATCH --export=ALL
#SBATCH --dependency=singleton

forcing_shape=$(parsed[:shape])
ncpu=$(parsed[:ncpu])
nthread=$(parsed[:nthread])
time_label=\$( date +"%Y%m%d%H%M%S" )

export SLURM_EXPORT_ENV=ALL
local_scratch=/scratch/\${USER}/job_\${SLURM_JOBID}


#export JULIA_NUM_THREADS=\$nthread
#echo "ncpu    = \$ncpu"
#echo "nthread = \$nthread"

echo "Local scratch    : \$local_scratch"
echo "SLURM_SUBMIT_DIR : \$SLURM_SUBMIT_DIR"

#export PATH="\$HOME/julia-1.6.3/bin:\$PATH"

if [ "\${SLURM_SUBMIT_DIR}" == "" ]; then
    wdir=`pwd`
else
    wdir=\${SLURM_SUBMIT_DIR}
fi


for i in \$(seq 1 10) ; do

    echo "The \${i}-th batch scan."

    julia \$wdir/../../../lib/continuation/01_main.jl \\
        --cont-varname gamma                         \\
        --method continuation                        \\
        --config \$wdir/config.jl                    \\
        --output-dir output_redo_\${forcing_shape}   \\
        --forcing-shape \$forcing_shape              \\
        --param $(parsed[:p_name])                   \\
        --param-value $(parsed[:p_val])              \\
        --scan-dir $scan_dir                         \\
        --scan-label "lb$(parsed[:ds_exp_lb])"       \\
        --param-factor $(parsed[:p_fac])             \\
        --ds-exp-lb $( - parsed[:ds_exp_lb] )        \\
        --scan-count 50


    if [ \$? -ne 0 ]; then 
        echo "Error occurs. Exit the loop."
        break
    fi 
done
"""
    )
end
