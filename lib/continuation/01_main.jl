using Distributed
using Formatting
using DataStructures
using JSON
using ArgParse

wdir = @__DIR__

println("Current directory: $wdir")

 
s = ArgParseSettings()
@add_arg_table s begin

    "--cont-varname"
        help = "The variable used for continuation. Valid input: `gamma`, `xi`."
        arg_type = String
        required = true

    "--method"
        help = "The scanning method. Available choices: `continuation`, `natural`."
        arg_type = String
        required = true

    "--output-dir"
        help = "Root output directory. It contains the initial snapshot file."
        arg_type = String
        required = true

    "--scan-label"
        help = "Output will be put into a subdirectory named with `--scan-label`."
        arg_type = String
        default = "noname"
        
    "--config"
        help = "Default configuration file."
        arg_type = String
        default = "config.jl"

    "--forcing-shape"
        help = "Freshwater forcing shape."
        arg_type = String
        required = true
 
    "--param"
        help = "Parameter that needs to vary. Currently accept cva_delta, Kh, Kv_iso, basin_width"
        arg_type = String
        required = true

    "--param-value"
        help = "Parameter value."
        arg_type = Float64

    "--param-factor"
        help = "This factor is going to multiply on the parameter so folder name can be pretty"
        arg_type = Float64
        default = 1e8


    # This section is for continuation scan
    "--scan-dir"
        help = "It accpets `pos`, `neg`."
        arg_type = String
        default = "pos"

    "--init-cond"
        help = "Initial condition snapshot file."
        arg_type = String
        default = ""

    "--ds-exp-lb"
        help = "."
        arg_type = Float64
        default = -5.0

    "--scan-count"
        help = "."
        arg_type = Int64
        default = 500

    "--newton-min-iter"
        help = "Newton method param."
        arg_type = Int64
        default = 3 

    "--newton-max-iter"
        help = "Newton method param."
        arg_type = Int64
        default = 50 

    "--newton-res-target"
        help = "Newton method param."
        arg_type = Float64
        default = 1e-13
    
    # This section is for natural scan
    "--scan-half-time"
        help = "The time required to scan Q from one bound to the other. So the total scan model time is double of this value. Unit: yr"
        arg_type = Int64
        default = 100
 
    "--record-interval"
        help = "The interval between two consecutive records. Unit: year."
        arg_type = Int64
        default = 10
 
    "--cva-off"
        help = "If set then cva is turned off by setting Kc = Kiso"
        action = :store_true
 

end
parsed = parse_args(ARGS, s)



println("Received parameters: ")
JSON.print(parsed, 4)

include(parsed["config"])

println("Configuration using: ")
JSON.print(cfg, 4)

casenames = Array{String}[]

root_output_dir = parsed["output-dir"]
mkpath(root_output_dir)

open(joinpath(root_output_dir, "config.json"), "w") do io
    JSON.print(io, cfg, 4)
end

function runCmds(cmds :: Array, filename::String = "/dev/stdout")
    for (i, cmd) in enumerate(cmds)
        if cmd == nothing
            println("The $(i)-th command is `nothing`. Skip this one.")
        else
            println(">> ", string(cmd))
            run(cmd) 
        end
    end
end

jobs = Channel(100)
completed_jobs = Channel(100)

if ! ( parsed["method"] in ["continuation", "natural"] )
    throw(ErrorException("Unknown method: $(parsed["method"])."))
end

Kh      = cfg["Kh"]
Kv_iso  = cfg["Kv_iso"]
cva_Δ   = cfg["cva_Δ"]
ξ       = cfg["ξ0"]
γ       = cfg["γ0"]
basin_width = cfg["domain_size"].basin
trans_width = cfg["Δϕ_trans"]
MLT_T = cfg["MLT_T"]
MLT_S = cfg["MLT_S"]
μ     = cfg["μ"]
spinup_snapshot_file = "spinup-snapshot.nc"
T_restore_days = cfg["T_restore_days"]

if parsed["param"] == "Kh"
    Kh      = parsed["param-value"]
elseif parsed["param"] == "Kv_iso"
    Kv_iso  = parsed["param-value"]
elseif parsed["param"] == "cva_delta"
    cva_Δ   = parsed["param-value"]
elseif parsed["param"] == "xi"
    ξ = parsed["param-value"]
elseif parsed["param"] == "gamma"
    γ = parsed["param-value"]
elseif parsed["param"] == "basin_width"
    basin_width = parsed["param-value"]
elseif parsed["param"] == "trans_width"
    trans_width = parsed["param-value"]
elseif parsed["param"] == "MLT_T"
    MLT_T = parsed["param-value"]
elseif parsed["param"] == "MLT_S"
    MLT_S = parsed["param-value"]
elseif parsed["param"] == "mu"
    μ = parsed["param-value"]
elseif parsed["param"] == "T_restore_days"
    T_restore_days = parsed["param-value"]
else
    throw(ErrorException("Unrecognized param: " * parsed["param"]))
end


function mkInitGuessCmd(
    output_dir      :: String,
    spinup_years    :: Int64;
    γ0  :: Union{Float64, Nothing} = nothing,
    ξ0  :: Union{Float64, Nothing} = nothing,
    transient_ξ     :: Float64 = 999.0
)

    if γ0 == nothing
        γ0 = γ
        println("Use γ0 = $γ0")
    end

    if ξ0 == nothing
        ξ0 = ξ
        println("Use ξ0 = $ξ0")
    end

    use_transient_ξ = transient_ξ != 999.0

    return `julia $wdir/02_gen_first_snapshot.jl
                --epsilon        $(cfg["ϵ"])
                --lat            $(cfg["lat"][1]) $(cfg["lat"][2])
                --Ny             $(cfg["Ny"])
                --Nz             $(cfg["Nz"])
                --MLT-T          $MLT_T
                --MLT-S          $MLT_S
                --MLT-shape      $(cfg["MLT_shape"])
                --Q-shape        $(parsed["forcing-shape"])
                --output-dir     $output_dir
                --snapshot-output $spinup_snapshot_file 
                --Q-forcing      $γ0
                --spinup-years   $spinup_years
                --days-per-record 360
                --T-restore-days $T_restore_days
                --basin-wide     $basin_width       
                --boundary-wide  $(cfg["domain_size"].bnd)       
                --cva            $(( parsed["cva-off"] ) ? "false" : "true" )
                --cva-softness   $cva_Δ
                --xi             $ξ0
                --Kh             $Kh
                --Kv-iso         $Kv_iso
                --init-cond      $(parsed["init-cond"])
                --trans-lat $(cfg["ϕc"])
                --trans-width $trans_width
                --use-transient-xi $use_transient_ξ
                --transient-xi $transient_ξ
                --transient-xi-duration $(cfg["transient_ξ_duration"])
                --mu             $μ
        `
end
        
casename_raw  = format("{:s}{:08d}", parsed["param"], Int64(round(parsed["param-value"]*parsed["param-factor"])))

cmds = []

if parsed["method"] == "continuation"
    
    if parsed["cont-varname"] == "gamma"
        varname = "γ"
        cont_file = "03_scan_cont_gamma.jl" 
    elseif parsed["cont-varname"] == "xi"
        varname = "ξ"
        cont_file = "11_scan_cont_xi.jl" 
    else
        throw(ErrorException("Unknown variable name: $(parsed["cont-varname"])"))
    end
   
    scan_dir = parsed["scan-dir"]
    scan_rng = cfg["$(varname)_rng_$(scan_dir)"]

    if scan_dir == "pos"
        scan_sign = 1.0
        begin_value = scan_rng[1]
    elseif scan_dir == "neg"
        scan_sign = -1.0
        begin_value = scan_rng[2]
    else
        throw(ErrorException("Unknown scan direction $scan_dir"))
    end

    casename  = format("CM_{:s}_{:s}", casename_raw, scan_dir)

    caseoutput = joinpath(root_output_dir, casename)
    mkpath(caseoutput)

    cmd_init_guess = nothing
    if isfile(joinpath(caseoutput, spinup_snapshot_file))
        println("File $spinup_snapshot_file already exists. Skip spinup.")
    else
        
        cmd_init_guess = mkInitGuessCmd(
            caseoutput,
            cfg["spinup"][scan_dir];
            Dict(
                Symbol("$(varname)0") => begin_value,
                :transient_ξ => ( ( cfg["use_transient_ξ_$(scan_dir)"] ) ? cfg["transient_ξ"] : 999.0 ) 
            )...
        )
    end

    cmd_cont = `julia $wdir/$cont_file
        --scan-dir  $scan_sign
        --scan-range $(scan_rng[1]) $(scan_rng[2])
        --scales $(cfg["scales"][1]) $(cfg["scales"][2]) $(cfg["scales"][3]) 
        --init-snapshot $(joinpath(caseoutput, spinup_snapshot_file))
        --casename $(parsed["scan-label"]) 
        --ds-exp-lb $(parsed["ds-exp-lb"])
        --scan-count $(parsed["scan-count"])
        --newton-min-iter $(parsed["newton-min-iter"])
        --newton-max-iter $(parsed["newton-max-iter"])
        --newton-res-target $(parsed["newton-res-target"])
    `

    push!(cmds, cmd_init_guess, cmd_cont)

elseif parsed["method"] == "natural"

    casename  = format("NA_{:s}", casename_raw)

    caseoutput = joinpath(root_output_dir, casename)
    mkpath(caseoutput)

    cmd_init_guess = nothing
    if isfile(joinpath(caseoutput, spinup_snapshot_file))
        println("File $spinup_snapshot_file already exists. Skip spinup.")
    else
        cmd_init_guess = mkInitGuessCmd(caseoutput, cfg["spinup-natural"], 0.0)
    end

    cmd_natural = `julia $(@__DIR__)/04_scan_natural.jl
        --scan-range $(cfg["F_rng_natural"][1]) $(cfg["F_rng_natural"][2])
        --scan-half-time $(parsed["scan-half-time"]) 
        --init-snapshot $(joinpath(caseoutput, spinup_snapshot_file))
        --casename $(parsed["scan-label"]) 
        --record-interval $(parsed["record-interval"])
        --num-records $(parsed["scan-count"])
    `
    push!(cmds, cmd_init_guess, cmd_natural)

end

exec_time = @elapsed runCmds(cmds)

println(format("It takes {:.2f} s to complete.", exec_time))
println("Program ends.")
