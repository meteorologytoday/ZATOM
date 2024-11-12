using Formatting
using ArgMacros

parsed = @dictarguments begin
    @argumentrequired Float64 γ          "--gamma"
    @argumentrequired Int64   ds_exp_lb   "--ds-exp-lb"
    @argumentrequired String  scan_dir    "--scan-dir"
    @argumentrequired String  shape       "--shape"
    @argumentrequired Int64   ncpu        "--ncpu"
    @argumentrequired Int64   nthread     "--nthread"
end

ds_exp_lb = abs(parsed[:ds_exp_lb])
scan_dir = parsed[:scan_dir]
shape = parsed[:shape]

if ! ( scan_dir in ["pos", "neg"] )
    throw(ErrorException("Unknown --scan-dir : $scan_dir"))
end

label = format("Xi60{:s}_P{:.3f}_lb{:d}_{:s}", parsed[:shape], parsed[:γ], ds_exp_lb, scan_dir)
casename=label

open("submit_$(label).sh", "w") do io
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

forcing_shape=$(shape)
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
        --cont-varname xi \\
        --method continuation                       \\
        --config \$wdir/config.jl     \\
        --output-dir output_redo_\${forcing_shape}       \\
        --forcing-shape \$forcing_shape             \\
        --param gamma                                  \\
        --param-value $(parsed[:γ])                  \\
        --scan-dir $scan_dir                   \\
        --scan-label "lb$(ds_exp_lb)"         \\
        --param-factor 1000                    \\
        --ds-exp-lb $(-ds_exp_lb)             \\
        --scan-count 50


    if [ \$? -ne 0 ]; then 
        echo "Error occurs. Exit the loop."
        break
    fi 
done
"""
    )
end
