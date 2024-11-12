using NCDatasets
using Formatting
using ArgParse
using JSON
using Statistics
using LinearAlgebra

include("helper_functions.jl")

"""
This function compute the Q of the following function.

Q(t) =     Q_rate * t                        if 0      <= t <  T_half
           Q_rate * (2 * T_half - t)         if T_half <= t <  2 * T_half
           0                                 otherwise

such that it is a triangular hump with the peak of (Q_rate * t) at T_half

"""
function Q_func(t, Q_rate, T_half)

    local ret_Q

    #println("Q_rate = $(Q_rate * 86400 * 360 / 1e6), T_half = $(T_half / 86400 / 360)")

    if 0 <= t < T_half
        ret_Q = Q_rate * t
    elseif t < 2 * T_half
        ret_Q = Q_rate * (2 * T_half - t)
    else
        ret_Q = 0.0
    end

    return ret_Q
end


println("""
This program takes a snapshot file (`--input-snapshot`) to begin continuition. The initial
direction is given by `--scan-dir`. Each scan contains `--scan-count` points. 

The output file are in the folder `[CASENAME]` and file name is "dddd.nc" CASENAME
is provided as `--casename`.

On initialization, if filename  "dddd.nc" is found, then by default the last record of
the file with the largest 4 digits number is used to initialize the model and continue
the scanning; if not, then start the entire new scan from "0001.nc"

""")

s = ArgParseSettings()
@add_arg_table s begin

    "--scan-range"
        help = "Scan range of Q in unit of Sv. Scan stops when Q is out of this range."
        nargs = 2
        arg_type = Float64
        required = true
 
    "--scan-half-time"
        help = "The time required to scan Q from one bound to the other. So the total scan model time is double of this value. Unit: yr"
        arg_type = Int64
        required = true
 
    "--record-interval"
        help = "The interval between two consecutive records. Unit: year."
        arg_type = Int64
        required = true
 
    "--num-records"
        help = "The number of records. This number is the dimension 'Ns' in the output file."
        arg_type = Int64
        required = true
 
    "--init-snapshot"
        help = "The snapshot file to begin with."
        arg_type = String
        required = true
 
    "--casename"
        help = "The casename of this scan. The output folder will be created using this name in the same folder as --input-snapshot"
        arg_type = String
        required = true
 
end


parsed = parse_args(ARGS, s)

JSON.print(parsed, 4)

sec_per_year = 360 * 86400.0

Q_min = parsed["scan-range"][1] * 1e6
Q_max = parsed["scan-range"][2] * 1e6

if Q_min > Q_max
    throw(ErrorException("Values `--scan-range` should be increasing."))
end

Q_rate = (Q_max - Q_min) / (sec_per_year * parsed["scan-half-time"])

println("(Q_min, Q_max) = ($Q_min, $Q_max)")
println("scan half time: $(parsed["scan-half-time"]) years.")
println(format("scan rate: {:.2e} Sv/yr", Q_rate / 1e6 * sec_per_year))


wrapped_Q_func = (t,) -> Q_func(t, Q_rate, parsed["scan-half-time"] * sec_per_year)


include("ZAOM/julia/constants.jl")
include("ZAOM/julia/ZAOM.jl")

using .ZAOM

function execPhysicalScan(;
    snapshot_file   :: String,         # This file is only used to construct the model
    init_file       :: String,         # This file will contain _X, Q, model_time
    output_file     :: String,
    Q_func          :: Function,
    t_beg           :: Float64,
    t_end           :: Float64,
    record_interval :: Float64,
    max_records     :: Int64,
)



    local model = ZAOM.loadModel(snapshot_file)
    local model_time
    local Q_now

    if t_beg > t_end
        throw(ErrorException("Error: t_beg > t_end"))
    end
#=
    function _callback(cmi)

        setXP!(cmi; s=cmi.s, x=ghost[:_X], p=ghost[:Q])
        model.state._X .= ghost[:_X]
        ZAOM.updateB!(model)
        ZAOM.updateSaltTransForcing!(model, ghost[:Q][1])
        ZAOM.diagnoseTransport!(model)
        
        setXP!(cmi; s=cmi.ṡ, x=ghost[:_dXds], p=ghost[:dQds])
    end
=#
    
    # Loading initial state
    Dataset(init_file, "r") do ds

        Ns = ds.dim["Ns"]        

        model.state._X .= reshape(ds["X_"][:, :, :, :, Ns], :)
        model_time = ds["model_time"][Ns]
        Q_now = Q_func(model_time)

        ZAOM.updateB!(model)
        ZAOM.updateSaltTransForcing!(model, Q_now)
       
    end

    if ! ( t_beg <= model_time <= t_end )
        println("Now model_time = $(model_time) is not within range [$t_beg, $t_end]")
        println("Not doing anything. Stop program.")
        return 1
    end
        
    total_steps = ceil(Int64, (t_end - model_time) / model.env.Δt )
    steps_per_record = floor(Int64, record_interval / model.env.Δt)

    println("Init_file: ", init_file)
    println(format("Begin time           : {:.2f}", t_beg))
    println(format("End   time           : {:.2f}", t_end))
    println(format("Model time           : {:.2f}", model_time))
    println(format("Model time step Δt   : {:.2f}", model.env.Δt))
    println(format("Record interval time :  {:.2f}", record_interval))
    println(format("Steps per record     : {:d}",   steps_per_record))
    println(format("Total steps needed to complete: {:d}", total_steps))
    println(format("Max records          : {:d}, equivalent to {:d} steps: ", max_records, max_records * steps_per_record))

    total_steps = min(total_steps, max_records * steps_per_record)
    println(format("Adjusted total_steps: {:d}", total_steps))


    println("Output file: $(output_file)")
    Dataset(output_file, "c") do ds

        local gd = model.env.gd_bib
        
        defDim(ds, "Ns", Inf)
        defDim(ds, "Nz", gd.Nz)
        defDim(ds, "Ny", gd.Ny)
        defDim(ds, "Nx",  2   )
        defDim(ds, "NX", model.env.NX)
        
        defDim(ds, "Nzp1", gd.Nz+1)
        defDim(ds, "Nyp1", gd.Ny+1)
        defDim(ds, "Nxp1", 3)

        ds.attrib["dlambda_b"] = model.env.Δλb
        ds.attrib["dlambda_i"] = model.env.Δλi

        # coordinate variables 
        z_T    = defVar(ds, "z_T", Float64, ("Nz",))
        z_W    = defVar(ds, "z_W", Float64, ("Nzp1",))
        y_T    = defVar(ds, "y_T", Float64, ("Ny",))
        y_V    = defVar(ds, "y_V", Float64, ("Nyp1",))
        Δx_T   = defVar(ds, "dx_T", Float64, ("Ny", "Nx",))

        z_T[:] = model.env.gd_bib.z_T[:, 1, 1] 
        z_W[:] = model.env.gd_bib.z_W[:, 1, 1]
        y_T[:] = model.env.gd_bib.ϕ_T[1, :, 1]
        y_V[:] = model.env.gd_bib.ϕ_V[1, :, 1]
        Δx_T[:] = model.env.gd_bib.Δx_T[1, :, :]

        #
        X_SS_   = defVar(ds, "X_SS_", Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))

        # simulation variables            

        Ψb      = defVar(ds, "Psib",   Float64, ("Nzp1", "Nyp1", "Ns"))
        X_      = defVar(ds, "X_",     Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
        bw      = defVar(ds, "bw",     Float64, ("Nz", "Ny", "Ns"))
        be      = defVar(ds, "be",     Float64, ("Nz", "Ny", "Ns"))
        diff_b  = defVar(ds, "diff_b", Float64, ("Nz", "Ny", "Ns"))
        ww      = defVar(ds, "ww",     Float64, ("Nzp1", "Ny", "Ns"))
        we      = defVar(ds, "we",     Float64, ("Nzp1", "Ny", "Ns"))
        vw      = defVar(ds, "vw",     Float64, ("Nz", "Nyp1", "Ns"))
        ui      = defVar(ds, "ui",     Float64, ("Nz", "Ny", "Ns"))
        Q       = defVar(ds, "Q",      Float64, ("Ns",))
        var_model_time = defVar(ds, "model_time", Float64, ("Ns",))
        
        # Diagnostic information
        ∫X      = defVar(ds, "int_X",    Float64, ("NX", "Ns",))
        ∫ADV    = defVar(ds, "int_ADV",  Float64, ("NX", "Ns",))
        ∫HDIF   = defVar(ds, "int_HDIF", Float64, ("NX", "Ns",))
        ∫VDIF   = defVar(ds, "int_VDIF", Float64, ("NX", "Ns",))
        ∫SS     = defVar(ds, "int_SS",   Float64, ("NX", "Ns",))
        ∫ERR    = defVar(ds, "int_ERR",  Float64, ("NX", "Ns",))
        ρc∫vT   = defVar(ds, "heat_transport", Float64, ("Nyp1", "Ns"))
        ρ∫vS    = defVar(ds, "salt_transport", Float64, ("Nyp1", "Ns"))

        ∫Δv = model.env.gd_bib.∫Δv
        function record!(_t)

            X_SS_[:, :, :, :, _t] = model.state.X_SS_
            Ψb[:, :, _t]      = model.state.Ψb
            X_[:, :, :, :, _t] = model.state.X_[:, :, :, :]
            bw[:, :, _t]      = model.state.b[:, :, 1]
            be[:, :, _t]      = model.state.b[:, :, 2]
            diff_b[:, :, _t]  = model.state.b[:, :, 2] - model.state.b[:, :, 1]
            ww[:, :, _t]      = model.state.w[:, :, 1]
            we[:, :, _t]      = model.state.w[:, :, 2]
            vw[:, :, _t]      = model.state.v[:, :, 1]
            ui[:, :, _t]      = model.state.u[:, :, 2]
            Q[_t]             = model.env.Q
            var_model_time[_t] = model_time

            ∫X[:, _t]         = model.state.diag[:∫X]     / ∫Δv
            ∫ADV[:, _t]       = model.state.diag[:∫ADV]   / ∫Δv
            ∫HDIF[:, _t]      = model.state.diag[:∫HDIF]  / ∫Δv
            ∫VDIF[:, _t]      = model.state.diag[:∫VDIF]  / ∫Δv
            ∫SS[:, _t]        = model.state.diag[:∫SS]    / ∫Δv
            ∫ERR[:, _t]       = model.state.diag[:∫ERR]   / ∫Δv
            
            ρc∫vT[:, _t]      = model.state.diag[:ρc∫vT]        
            ρ∫vS[:, _t]       = model.state.diag[:ρ∫vS]        
            
            println(format("Record the {:d}-th result. ", _t))
        end

        # Compute time steps

            
        println("Starting stepping the model...")
        for s = 1:total_steps

            Q_now = Q_func(model_time)

            print(format("\rQ_now: {:f} Sv ({:d}/{:d})", Q_now / 1e6, s, total_steps))

            ZAOM.updateB!(model)
            ZAOM.updateSaltTransForcing!(model, Q_now)
            ZAOM.diagnoseTransport!(model)
                
            ZAOM.stepModel_rk2!(model)
            model_time += model.env.Δt

            if mod(s, steps_per_record) == 0
                println()
                println(format("Current model time: {:.2f}, Q_now = {:.2f} Sv.  ({:d}/{:d})", model_time, Q_now / 1e6, s, total_steps))
                record!(Int64(s/steps_per_record))
            end

        end
    end

    println("Function `execPhysicalScan` Done.")

    return 0

end

println("Start Program")
println("Detecting initialization...")

while true

    global status, init_file, next_file = detectInitialization(
        dirname(parsed["init-snapshot"]),
        parsed["casename"],
    )

    global next_file_tmp = joinpath(dirname(next_file), "$(basename(next_file)).tmp")

    if status == :INIT

        println("No previous scan. Create a new one: $(next_file_tmp)")

        global model = ZAOM.loadModel(parsed["init-snapshot"])
        local model_time = 0.0
        local Q0 = wrapped_Q_func(model_time)

        println("Q0 = ", Q0/1e6, " Sv")
        ZAOM.updateSaltTransForcing!(model, Q0)

        println("Outputting: $(next_file_tmp)") 
        Dataset(next_file_tmp, "c") do ds

            local gd = model.env.gd_bib
            
            defDim(ds, "Ns", Inf)
            defDim(ds, "Nz", gd.Nz)
            defDim(ds, "Ny", gd.Ny)
            defDim(ds, "Nx",  2   )
            defDim(ds, "NX", model.env.NX)
 
            X_      = defVar(ds, "X_",      Float64, ("Nz", "Ny", "Nx", "NX", "Ns"))
            Q       = defVar(ds, "Q",       Float64, ("Ns",))
            var_model_time = defVar(ds, "model_time",       Float64, ("Ns",))
           
            X_[:, :, :, :, 1]      = model.state.X_
            Q[1]                   = Q0
            var_model_time[1]      = model_time
        end

        mv(next_file_tmp, next_file; force=true)

    elseif status == :CONT
       
        println("Initial file determined") 
        break

    end

end

println("Init the scan with : $init_file")
println("Running natural scanning...")

status = execPhysicalScan(;
    snapshot_file     = parsed["init-snapshot"],
    init_file         = init_file,
    output_file       = next_file_tmp,
    Q_func            = wrapped_Q_func,
    t_beg             = 0.0,
    t_end             = 2 * parsed["scan-half-time"] * sec_per_year,
    record_interval   = parsed["record-interval"] * sec_per_year,
    max_records       = parsed["num-records"],
)

if status == 0
    println("Moving $(next_file_tmp) to $(next_file)")        
    mv(next_file_tmp, next_file; force=true)
else
    println("Error: execPhysicalScan returns nonzero code $(status). Please check.")
end


println("End Program.")




