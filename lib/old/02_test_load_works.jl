using NCDatasets
using Formatting
include("ZAOM/constants.jl")
include("ZAOM/ZAOM.jl")

rec_output_file = "output03_record_testloading.nc"

model = ZAOM.loadModel("output02_snapshot.nc")

ZAOM.init!(model)


@time Dataset(rec_output_file, "c") do ds

    gd = model.env.gd_bib

    defDim(ds, "Nz", gd.Nz)
    defDim(ds, "Ny", gd.Ny)
    defDim(ds, "Nx",   2)
    defDim(ds, "NX", model.env.NX)
    
    defDim(ds, "Nzp1", gd.Nz+1)
    defDim(ds, "Nyp1", gd.Ny+1)
    defDim(ds, "Nxp1", 3)
    
    defDim(ds, "time", Inf)
   
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
    X_SS_   = defVar(ds, "X_SS_", Float64, ("Nz", "Ny", "Nx", "NX",))
    X_SS_[:, :, :, :] = model.state.X_SS_
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
    
    # Diagnostic information
    ∫X      = defVar(ds, "int_X",    Float64, ("NX", "time",))
    ∫ADV    = defVar(ds, "int_ADV",  Float64, ("NX", "time",))
    ∫HDIF   = defVar(ds, "int_HDIF", Float64, ("NX", "time",))
    ∫VDIF   = defVar(ds, "int_VDIF", Float64, ("NX", "time",))
    ∫SS     = defVar(ds, "int_SS",   Float64, ("NX", "time",))
    ∫ERR    = defVar(ds, "int_ERR",  Float64, ("NX", "time",))


    ∫Δv = model.env.gd_bib.∫Δv
    function record!(_t)

        Ψb[:, :, _t]      = model.state.Ψb
        X_[:, :, :, :, _t] = model.state.X_[:, :, :, :]
        bw[:, :, _t]      = model.state.b[:, :, 1]
        be[:, :, _t]      = model.state.b[:, :, 2]
        diff_b[:, :, _t]  = model.state.b[:, :, 2] - model.state.b[:, :, 1]
        ww[:, :, _t]      = model.state.w[:, :, 1]
        we[:, :, _t]      = model.state.w[:, :, 2]
        vw[:, :, _t]      = model.state.v[:, :, 1]
        ui[:, :, _t]      = model.state.u[:, :, 2]
        ∫X[:, _t]         = model.state.diag[:∫X]     / ∫Δv
        ∫ADV[:, _t]       = model.state.diag[:∫ADV]   / ∫Δv
        ∫HDIF[:, _t]      = model.state.diag[:∫HDIF]  / ∫Δv
        ∫VDIF[:, _t]      = model.state.diag[:∫VDIF]  / ∫Δv
        ∫SS[:, _t]        = model.state.diag[:∫SS]    / ∫Δv
        ∫ERR[:, _t]       = model.state.diag[:∫ERR]   / ∫Δv
        
    end

    for t=1:360*2

        print(format("Step: {:d}\n", t))
        ZAOM.stepModel_rk2!(model)

        record!(t)

    end

end

println("End Program.")




