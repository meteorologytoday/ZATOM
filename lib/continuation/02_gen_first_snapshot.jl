using NCDatasets
using Formatting
using ArgParse
using JSON

include("../../julia/constants.jl")
include("../../julia/ZAOM.jl")

s = ArgParseSettings()
@add_arg_table s begin

    "--output-dir"
        help = "Output dir"
        arg_type = String
        required = true

    "--Q-forcing"
        help = "Freshwater flux forcing in Sv."
        arg_type = Float64
        required = true

    "--Q-shape"
        help = "Shape of Q. Currerntly available: flat, cosine."
        arg_type = String
        required = true


    "--xi"
        help = "Value of ξ. "
        arg_type = Float64
        required = true
 
    "--use-transient-xi"
        help = "Use transient ξ to help spin up (--transient-xi and --transient-xi-duration). "
        arg_type = Bool
        default = false
 
    "--transient-xi"
        help = "Init value of ξ. "
        arg_type = Float64
        default = NaN
 
    "--transient-xi-duration"
        help = "Duration of the --transient-xi in years."
        arg_type = Int64
        default = 10
 
    "--cva"
        help = ""
        arg_type = Bool
        default = true
 
    "--cva-softness"
        help = ""
        arg_type = Float64
        required = true
 
    "--T-restore-days"
        help = "Restoring timescale of temperature (days)."
        arg_type = Float64
        default  = 5.0


    "--spinup-years"
        help = "Spinup time in years."
        arg_type = Int64
        default  = 5

    "--days-per-record"
        help = "In days."
        arg_type = Int64
        default = 30
 
    "--basin-wide" 
        help = "Basin wide in degrees."
        arg_type = Float64
        default  = 20.0

    "--boundary-wide" 
        help = "Boundary wide in degrees."
        arg_type = Float64
        default  = 4.0

    "--record-output" 
        help = "Record output filename"
        arg_type = String
        default  = "output-spinup-record.nc"

    "--snapshot-output" 
        help = "Snapshot output filename"
        arg_type = String
        default  = "output-spinup-snapshot.nc"

    "--Kh"
        help = "Horizontal diffusivity."
        arg_type = Float64
        default  = 4.4e4

    "--Kv-iso"
        help = "Vertical diffusivity."
        arg_type = Float64
        default  = 1e-4

    "--epsilon"
        help = "Rayleigh friction."
        arg_type = Float64
        required = true

    "--Ny"
        help = "Number of points in y direction."
        arg_type = Int64
        default  = 30

    "--Nz"
        help = "Number of points in z direction."
        arg_type = Int64
        default  = 30

    "--lat"
        help = "North and south boundary in latitude."
        arg_type = Float64
        nargs = 2
        required = true

    "--MLT-T"
        help = "Mixed-layer thickness for temperautre."
        arg_type = Float64
        default  = 100.0

    "--MLT-S"
        help = "Mixed-layer thickness for salinity."
        arg_type = Float64
        default  = 100.0


    "--MLT-shape"
        help = "Mixed-layer function shape. It can be `linear` or `step`."
        arg_type = String
        default  = "linear"


    "--init-cond"
        help = "The snapshot file used as initial condition. Only temperature and salinity are used."
        arg_type = String
        default = ""

    "--trans-lat"
        help = "Transition latitude (in degrees) of fresh water forcing."
        arg_type = Float64
        required = true

    "--trans-width"
        help = "Trnsition width (in degress) of fresh water forcing."
        arg_type = Float64
        required = true

    "--cont-phy-years"
        help = "The physically stepping years after steady state is found by newton method."
        arg_type = Int64
        default = 10

    "--dT-east"
        help = "The temperature difference along the eastern boundary."
        arg_type = Float64
        default = 25.0

    "--dT-west"
        help = "The temperature difference along the western boundary."
        arg_type = Float64
        default = 25.0

    "--T-south"
        help = "The temperature along the southern boundary."
        arg_type = Float64
        default = 30.0

    "--mu"
        help = "The sensitivity of geostrophic balance."
        arg_type = Float64
        default = 1.0

end

parsed = parse_args(ARGS, s)



JSON.json(parsed; pretty=true)

# ===== [BEGIN setting] =====
days_per_year = 360

sim_days      =  days_per_year * parsed["spinup-years"] 
steps_per_day =    1
days_per_rec  =    parsed["days-per-record"]

do_newton = true

newton_steps_min = 5
newton_steps_max = 100
newton_tol = 1e-15

do_contphys = true
contphys_days = parsed["cont-phy-years"] * days_per_year

Δλ_deg  = parsed["basin-wide"]
Δλb_deg = parsed["boundary-wide"]


transient_ξ = parsed["xi"]
target_ξ   = parsed["xi"]
timestep_change_ξ = Int(0)

if parsed["use-transient-xi"]
    transient_ξ = parsed["transient-xi"]
    if parsed["transient-xi-duration"] < 0
        throw(ErrorException("--transient-xi-duration should be a positive integer. Now it is $(parsed["transient-xi-duration"])"))
    end
    timestep_change_ξ = Int(steps_per_day * parsed["transient-xi-duration"] * days_per_year) 
end

println("timestep_change_ξ = $(timestep_change_ξ)")

# ===== [END setting] =====

mkpath(parsed["output-dir"])
for varname in ["record-output", "snapshot-output"]
    parsed[varname] = joinpath(parsed["output-dir"], parsed[varname])
    println(format("parsed[{:s}] = \"{:s}\"", varname, parsed[varname]))
end

newton_steps_min *= do_newton
newton_steps_max *= do_newton

phys_steps = Int(sim_days      * steps_per_day)
continue_steps = Int(contphys_days * steps_per_day)


total_steps_max = phys_steps + newton_steps_max + do_contphys * continue_steps

if newton_steps_max < newton_steps_min 
    throw(ErrorException("newton_steps_max should be greater or equal to newton_steps_min"))
end

ϕn = deg2rad(parsed["lat"][2])
ϕs = deg2rad(parsed["lat"][1])

if ϕn <= ϕs
    throw(ErrorException("Option --lat should be lat_s and lat_n."))
end

env = ZAOM.Env(;
    ϵ  = parsed["epsilon"],
    Ω  = 7.3e-5,
    Ny = parsed["Ny"], #60,
    Nz = parsed["Nz"],# 128,
    Δλb = deg2rad(Δλb_deg),
    Δλi = deg2rad(Δλ_deg-Δλb_deg*2),
    ϕn = deg2rad(parsed["lat"][2]),
    ϕs = deg2rad(parsed["lat"][1]),
    H  = 4500.0,
    R  = 6400e3,
    Δt     = 86400.0 / steps_per_day,
    Kh     = parsed["Kh"],
    Kv     = parsed["Kv-iso"],
    Kv_cva = min(( (parsed["cva"]) ? 1e4 : 1.0 ) * parsed["Kv-iso"], 1.0),
    cva_Δ  = parsed["cva-softness"],
    tracer_forcings = [
        ZAOM.TracerForcing(:SURFACE, 86400.0 * parsed["T-restore-days"]),
        ZAOM.TracerForcing(:EVERYWHERE, 86400.0 * 360 * 1e6),
    ],
    ϕc = deg2rad(parsed["trans-lat"]),
    Δϕ_trans = deg2rad(parsed["trans-width"]),
    Q = parsed["Q-forcing"] * 1e6,
    Q_shape  = parsed["Q-shape"],
    ξ        = parsed["xi"],
    MLT_T    = parsed["MLT-T"],
    MLT_S    = parsed["MLT-S"],
    MLT_shape= parsed["MLT-shape"],
    μ        = parsed["mu"],
)

#env.Khx .= (5000 * exp.(env.gd_b.z_T ./ 100) .+ 0.0)[:, :, 1]


model = ZAOM.Model(env)
gd_bb = model.env.gd_bb

# initial condition
ΔT = .05 / (ZAOM.α*ZAOM.g)
ΔS = .05 / (ZAOM.β*ZAOM.g) * 0

ΔT_east = parsed["dT-east"]
ΔT_west = parsed["dT-west"]

S_init = similar(gd_bb.ϕ_T) 
S_init .= ZAOM.S0

SST = similar(gd_bb.ϕ_T)
SST[:, :, 1] .= ( parsed["T-south"] .+ parsed["dT-west"] / 2.0 * ( cos.( π * ( gd_bb.ϕ_T .- deg2rad(10.0) ) / deg2rad(60.0) ) .- 1 ))[:, :, 1] 
SST[:, :, 2] .= ( parsed["T-south"] .+ parsed["dT-east"] / 2.0 * ( cos.( π * ( gd_bb.ϕ_T .- deg2rad(10.0) ) / deg2rad(60.0) ) .- 1 ))[:, :, 2]


SST = SST[1:1, :, :]  # only need one layer

T_init = similar(gd_bb.ϕ_T)
T_init .= SST .* exp.(gd_bb.z_T ./ 450)

model.state.T .= T_init
model.state.S .= S_init

if parsed["init-cond"] != ""
    println("Using initial condition with file: $(parsed["init-cond"])")
    tmp_model = ZAOM.loadModel(parsed["init-cond"])
    model.state._X .= tmp_model.state._X 
    tmp_model = nothing
end


# Notice that target has to be z-independent
model.state.T_tgt[:, :, :] .= SST


# Initiate model (update buoyancy and such)
ZAOM.updateB!(model)
ZAOM.updateSaltTransForcing!(model)

∫Δv = model.env.gd_bib.∫Δv

steps_per_rec = steps_per_day * days_per_rec
first_time_newton = true
@time Dataset(parsed["record-output"], "c") do ds

    defDim(ds, "Nz", gd_bb.Nz)
    defDim(ds, "Ny", gd_bb.Ny)
    defDim(ds, "Nx",   2)
    defDim(ds, "NX", model.env.NX)
    
    defDim(ds, "Nzp1", gd_bb.Nz+1)
    defDim(ds, "Nyp1", gd_bb.Ny+1)
    defDim(ds, "Nxp1", 3)
    
    defDim(ds, "time", Inf)
    defDim(ds, "newton", newton_steps_max)
   
    #
    phase_idx = defVar(ds, "phase_idx", Float64, ("time",))

    # coordinate variables 
    z_T    = defVar(ds, "z_T", Float64, ("Nz",))
    z_W    = defVar(ds, "z_W", Float64, ("Nzp1",))
    y_T    = defVar(ds, "y_T", Float64, ("Ny",))
    y_V    = defVar(ds, "y_V", Float64, ("Nyp1",))
    Δx_T   = defVar(ds, "dx_T", Float64, ("Nx",))

    z_T[:] = model.env.gd_bib.z_T[:, 1, 1] 
    z_W[:] = model.env.gd_bib.z_W[:, 1, 1]
    y_T[:] = model.env.gd_bib.ϕ_T[1, :, 1]
    y_V[:] = model.env.gd_bib.ϕ_V[1, :, 1]
    Δx_T[:] = model.env.gd_bib.Δx_T[1, 1, :]

    #
    X_SS_        = defVar(ds, "X_SS_",       Float64, ("Nz", "Ny", "Nx", "NX", "time"))
    X_VDIFU_     = defVar(ds, "X_VDIFU_",    Float64, ("Nz", "Ny", "Nx", "NX", "time"))
    X_HDIFU_     = defVar(ds, "X_HDIFU_",    Float64, ("Nz", "Ny", "Nx", "NX", "time" ))
    X_ZNLHDIFU_  = defVar(ds, "X_ZNLHDIFU_", Float64, ("Nz", "Ny", "Nx", "NX", "time" ))
    X_CVA_       = defVar(ds, "X_CVA_",      Float64, ("Nz", "Ny", "Nx", "NX", "time" ))
    X_ADV_ZOC_   = defVar(ds, "X_ADV_ZOC_",  Float64, ("Nz", "Ny", "Nx", "NX", "time" ))
    X_ADV_MOC_   = defVar(ds, "X_ADV_MOC_",  Float64, ("Nz", "Ny", "Nx", "NX", "time" ))
    X_FRC_       = defVar(ds, "X_FRC_",      Float64, ("Nz", "Ny", "Nx", "NX", "time" ))
    dXdt_SUM_    = defVar(ds, "dXdt_SUM_",   Float64, ("Nz", "Ny", "Nx", "NX", "time" ))

    u_ZOC   = defVar(ds, "u_ZOC",     Float64, ("Nz", "Ny", "Nxp1", "time"))
    v_ZOC   = defVar(ds, "v_ZOC",     Float64, ("Nz", "Nyp1", "Nx", "time"))
    w_ZOC   = defVar(ds, "w_ZOC",     Float64, ("Nzp1", "Ny", "Nx", "time"))

    u_MOC   = defVar(ds, "u_MOC",     Float64, ("Nz", "Ny", "Nxp1", "time"))
    v_MOC   = defVar(ds, "v_MOC",     Float64, ("Nz", "Nyp1", "Nx", "time"))
    w_MOC   = defVar(ds, "w_MOC",     Float64, ("Nzp1", "Ny", "Nx", "time"))



    MLT_T_w   = defVar(ds, "MLT_T_w", Float64, ("Nz", "Ny", "Nx",))
    MLT_S_w   = defVar(ds, "MLT_S_w", Float64, ("Nz", "Ny", "Nx",))
    MLT_T_w[:, :, :] = ZAOM.genMixedlayerWeight(model.env.gd_bib.z_T, model.env.MLT_T, model.env.MLT_shape)
    MLT_S_w[:, :, :] = ZAOM.genMixedlayerWeight(model.env.gd_bib.z_T, model.env.MLT_S, model.env.MLT_shape)
 
    Khx = defVar(ds, "Khx", Float64, ("Nz", "Ny", ))
    Khx[:, :] = model.env.Khx
 
    # simulation variables            

    Ψb      = defVar(ds, "Psib",   Float64, ("Nzp1", "Nyp1", "time"))
    X_      = defVar(ds, "X_",     Float64, ("Nz", "Ny", "Nx", "NX", "time"))
    bw      = defVar(ds, "bw",     Float64, ("Nz", "Ny", "time"))
    be      = defVar(ds, "be",     Float64, ("Nz", "Ny", "time"))
    diff_b  = defVar(ds, "diff_b", Float64, ("Nz", "Ny", "time"))
    ww      = defVar(ds, "ww",     Float64, ("Nzp1", "Ny", "time"))
    we      = defVar(ds, "we",     Float64, ("Nzp1", "Ny", "time"))
    vw      = defVar(ds, "vw",     Float64, ("Nz", "Nyp1", "time"))
    ui      = defVar(ds, "ui",     Float64, ("Nz", "Ny", "time"))
    χ       = defVar(ds, "chi",    Float64, ("Nzp1", "Ny", "time"))
    ξ       = defVar(ds, "xi",     Float64, ("time",))
    qe      = defVar(ds, "qe",     Float64, ("Nz", "Ny", "NX", "time"))
    qw      = defVar(ds, "qw",     Float64, ("Nz", "Ny", "NX", "time"))


    
    
    residue = defVar(ds, "residue", Float64, ("newton",))
    
    # Diagnostic information
    ∫X      = defVar(ds, "int_X",    Float64, ("NX", "time",))
    ∫ADV    = defVar(ds, "int_ADV",  Float64, ("NX", "time",))
    ∫HDIF   = defVar(ds, "int_HDIF", Float64, ("NX", "time",))
    ∫VDIF   = defVar(ds, "int_VDIF", Float64, ("NX", "time",))
    ∫SS     = defVar(ds, "int_SS",   Float64, ("NX", "time",))
    ∫ERR    = defVar(ds, "int_ERR",  Float64, ("NX", "time",))
    ρc∫vT   = defVar(ds, "heat_transport", Float64, ("Nyp1", "time"))
    ρ∫vS    = defVar(ds, "salt_transport", Float64, ("Nyp1", "time"))

    function record!(_t, _phase)

        phase_idx[_t]     = Dict(
            :init         => 0.0,
            :phys     => 1.0,
            :newton       => 2.0,
            :contphys => 3.0,
            :finalize     => 4.0,
        )[_phase]


        X_SS_[:, :, :, :, _t]    = model.state.X_SS_
        X_VDIFU_[:, :, :, :, _t] = model.state.X_VDIFU_
        X_HDIFU_[:, :, :, :, _t] = model.state.X_HDIFU_
        X_ZNLHDIFU_[:, :, :, :, _t] = model.state.X_ZNLHDIFU_
        X_CVA_[:, :, :, :, _t] = model.state.X_CVA_
        X_ADV_ZOC_[:, :, :, :, _t] = model.state.X_ADV_ZOC_
        X_ADV_MOC_[:, :, :, :, _t] = model.state.X_ADV_MOC_
        X_FRC_[:, :, :, :, _t] = model.state.X_FRC_
        dXdt_SUM_[:, :, :, :, _t] = model.state.dXdt_SUM_

        u_ZOC[:, :, :, _t]  = model.state.u_ZOC
        v_ZOC[:, :, :, _t]  = model.state.v_ZOC
        w_ZOC[:, :, :, _t]  = model.state.w_ZOC
        u_MOC[:, :, :, _t]  = model.state.u_MOC
        v_MOC[:, :, :, _t]  = model.state.v_MOC
        w_MOC[:, :, :, _t]  = model.state.w_MOC
 

        Ψb[:, :, _t]      = model.state.Ψb
        X_[:, :, :, :, _t] = model.state.X_[:, :, :, :]
        bw[:, :, _t]      = model.state.b[:, :, 1]
        be[:, :, _t]      = model.state.b[:, :, 2]
        diff_b[:, :, _t]  = model.state.b[:, :, 2] - model.state.b[:, :, 1]
        ww[:, :, _t]      = model.state.w[:, :, 1]
        we[:, :, _t]      = model.state.w[:, :, 2]
        vw[:, :, _t]      = model.state.v[:, :, 1]
        ui[:, :, _t]      = model.state.u[:, :, 2]
        χ[:, :, _t]       = model.state.χ[:, :]
        ξ[_t]             = model.env.ξ
        qw[:, :, :, _t]   = model.state.diag[:q_cva][:, :, 1, :]
        qe[:, :, :, _t]   = model.state.diag[:q_cva][:, :, 2, :]
        ∫X[:, _t]         = model.state.diag[:∫X]     / ∫Δv
        ∫ADV[:, _t]       = model.state.diag[:∫ADV]   / ∫Δv
        ∫HDIF[:, _t]      = model.state.diag[:∫HDIF]  / ∫Δv
        ∫VDIF[:, _t]      = model.state.diag[:∫VDIF]  / ∫Δv
        ∫SS[:, _t]        = model.state.diag[:∫SS]    / ∫Δv
        ∫ERR[:, _t]       = model.state.diag[:∫ERR]   / ∫Δv

        ρc∫vT[:, _t]      = model.state.diag[:ρc∫vT]        
        ρ∫vS[:, _t]       = model.state.diag[:ρ∫vS]        
#        dbdy[:, :, _t]    = reshape(model.core.hds.sop.V_m∂y_T * model.state.B
    end

    ds.attrib["total_fwf"] = ZAOM.calTotalFreshwaterFlux(model)


    t = 1
    r_prev = 1.0

    phys_steps_in_use = 0
    phase = :init
    step = 0

    while true

        print(format("[{:s}] Step: {:d}/{:d} ({:.2f}%)\n", String(phase), step, total_steps_max, 100*step/total_steps_max))

        record_flag = false

        if phase == :init

            phase = :phys
            step = 0
            record_flag = true
            phys_steps_in_use = phys_steps
            
            ZAOM.updateSaltTransForcing!(model; ξ = transient_ξ)
            
        elseif phase == :phys || phase == :contphys

            if step < phys_steps_in_use

                if mod(step, steps_per_rec) == 0
                    record_flag = true
                end

                if phase == :phys && step == timestep_change_ξ 
                    ZAOM.updateSaltTransForcing!(model; ξ = target_ξ)
                    println("Transient time up. Update ξ to $(model.env.ξ).") 
                end


                ZAOM.stepModel_rk2!(model)
                #ZAOM.stepModel_newton!(model)

            else

                if phase == :phys
                    if do_newton
                        phase = :newton
                    elseif do_contphys
                        phase = :contphys
                    else
                        phase = :finalize
                    end
                elseif phase == :contphys
                    phase = :finalize
                else
                    throw(ErrorException("Unknown phase: " * string(phase)))
                end

                step = 0

            end

        elseif phase == :newton
        
            record_flag = true

            if first_time_newton 
                ZAOM.updateSaltTransForcing!(model; ξ = target_ξ)
                println("Update ξ to $(model.env.ξ)")
                global first_time_newton = false
            end 

            r, logDetJ = ZAOM.stepModel_newton!(model; verbose=true)
            #_, _, r = ZAOM.stepModel_newton!(model, res_target=1e-10, min_iter=5, max_iter=20, verbose=true, fastforward=true, )
            
            residue[step] = r

            println(format("Newton step {:d} , log(det(J)) = {:.2e}, sign = {:d}. |Δ| = {:.5e}. Rel chg: {:.2f}%.", step, logDetJ..., r, (r - r_prev) / r_prev * 100.0))
            
            if r > 1e5
                throw(ErrorException("Newton method explodes. Residue = $r"))
            end
 
            if step >= newton_steps_min && r < newton_tol 
                
                println("Converge criteria met. Break the loop...")

                phase = :contphys
                phys_steps_in_use = continue_steps
                step = 0 

            elseif newton_steps_min != 0 && step == newton_steps_max

                println("Maximum newton steps reached without convergence.")
                phase = :finalize 
                step = 0                

            end

            r_prev = r

        elseif phase == :finalize
            break
        end

        if record_flag
            record!(t, phase)
            t+=1
        end

        step += 1

    end

end

println("Saving snapshot...")
ZAOM.saveModel(model, parsed["snapshot-output"])
println("End Program.")
