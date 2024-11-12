function getCoordinateVariable(m::Model)
    return Dict(
#        "f"               => ( m.env.gi.c_f,                                                ("Nx", "Ny",) ),
    )
end

function getCompleteVariableList(m::Model)
    tmd = m.tmd_engine.state
    dyn = m.dyn_engine.state

    return Dict(
        "T"               => ( PermutedDimsArray(view(tmd.X, :, :, :, 1), (2, 3, 1)), ("Nx", "Ny", "Nz_f") ),
        "S"               => ( PermutedDimsArray(view(tmd.X, :, :, :, 2), (2, 3, 1)), ("Nx", "Ny", "Nz_f") ),
        "T_ML"            => ( view(tmd.X_ML, :, :, 1),          ("Nx", "Ny",) ),
        "S_ML"            => ( view(tmd.X_ML, :, :, 2),          ("Nx", "Ny",) ),
        "h_ML"            => ( tmd.h_ML,                            ("Nx", "Ny",) ),
        "swflx"           => ( tmd.swflx,                           ("Nx", "Ny",) ),
        "nswflx"          => ( tmd.nswflx,                          ("Nx", "Ny",) ),
        "Phi"             => ( dyn.Φ,                               ("Nx", "Ny",) ),
        "u_total"         => ( dyn.u_total,                         ("Nx", "Ny",   "Nz_c") ),
        "v_total"         => ( dyn.v_total,                         ("Nx", "Nyp1", "Nz_c") ),
    )
end


function getBasicRecorder(
    m :: Model,
)
    

    # Setting up recorder
    complete_varlist = getCompleteVariableList(m)
    varlist = []
    for varname in keys(complete_varlist)
        println(format("Using varaible: {:s}", varname))
        push!(varlist, ( varname, complete_varlist[varname]... ) )
    end

    coord_varlist = getCoordinateVariable(m)
    cvarlist = []
    for varname in keys(coord_varlist)
        println(format("Using varaible: {:s}", varname))
        push!(cvarlist, ( varname, coord_varlist[varname]... ) )
    end

    return RecordTool.Recorder(
        Dict(
            "NX"      => m.env.NX,
            "Nyp1"    => m.env.Ny+1,
            "Nz_fp1"  => m.env.Nz_f+1,
            "Nx"      => m.env.Nx,
            "Ny"      => m.env.Ny,
            "Nz_f"    => m.env.Nz_f,
            "Nz_c"    => m.env.Nz_c,
            "z_bnd_f" => length(m.env.z_bnd_f),
        ), varlist, Dict(),
        other_vars = cvarlist
    )
     

end

