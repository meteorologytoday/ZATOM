function saveModel(
    m :: Model,
    filename :: String,
)
    @fast_extract m

    gd = m.env.gd_bib

    Dataset(filename, "c") do ds

        defDim(ds, "Nz", gd.Nz)
        defDim(ds, "Ny", gd.Ny)
        defDim(ds, "Nx",   2)
        defDim(ds, "NX", ev.NX)
    
        defDim(ds, "Nzp1", gd.Nz+1)
        defDim(ds, "Nyp1", gd.Ny+1)
        defDim(ds, "Nxp1", 3)
 
        ds.attrib["Omega"]   = gd.Ω
        ds.attrib["epsilon"] = gd.ϵ
 
        ds.attrib["dlambda_b"] = ev.Δλb
        ds.attrib["dlambda_i"] = ev.Δλi
        
        ds.attrib["phi_n"] = gd.ϕn
        ds.attrib["phi_s"] = gd.ϕs
        
        ds.attrib["H"] = gd.H
        ds.attrib["R"] = gd.R
        
        ds.attrib["dt"]     = ev.Δt
        ds.attrib["Kh"]     = ev.Kh
        ds.attrib["Kv"]     = ev.Kv
        ds.attrib["Kv_cva"] = ev.Kv_cva
        ds.attrib["cva_delta"] = ev.cva_Δ
        
        ds.attrib["phi_c"]       = ev.ϕc
        ds.attrib["dphi_trans"]  = ev.Δϕ_trans
        ds.attrib["Q"]           = ev.Q
        ds.attrib["Q_shape"]     = ev.Q_shape
        ds.attrib["xi"]          = ev.ξ
        ds.attrib["MLT_T"]    = ev.MLT_T
        ds.attrib["MLT_S"]    = ev.MLT_S
        ds.attrib["MLT_shape"]    = ev.MLT_shape
        ds.attrib["mu"]    = ev.μ
        
        X_     = defVar(ds, "X_",      Float64, ("Nz", "Ny", "Nx", "NX"))
        X_tgt_ = defVar(ds, "X_tgt_",  Float64, ("Nz", "Ny", "Nx", "NX"))

        forcing_style = defVar(ds, "forcing_style",  Float64, ("NX",))
        forcing_τ     = defVar(ds, "forcing_time",   Float64, ("NX",))

        X_[:,:,:,:]     = st.X_
        X_tgt_[:,:,:,:] = st.X_tgt_

        fs = zeros(Float64, ev.NX) .- 1
        fτ = zeros(Float64, ev.NX) .- 1

        fs_mapping = Dict(
            :NONE       => 0.0,
            :SURFACE    => 1.0,
            :EVERYWHERE => 2.0,
        )

        for k=1:ev.NX
            tracer_forcing = ev.tracer_forcings[k]
            fs[k] = fs_mapping[tracer_forcing.style]
            fτ[k] = tracer_forcing.τ
        end
        forcing_style[:] = fs
        forcing_τ[:]     = fτ
    end
 
end

function loadModel(
    filename :: String,
)

    m = nothing

    fs_mapping = Dict(
        0.0 => :NONE,
        1.0 => :SURFACE,
        2.0 => :EVERYWHERE,
    )


    Dataset(filename, "r") do ds

        tracer_forcings = Array{TracerForcing}(undef, ds.dim["NX"])
        for k=1:length(tracer_forcings)
            tracer_forcings[k] = TracerForcing(fs_mapping[ds["forcing_style"][k]], ds["forcing_time"][k])
        end

        env =  Env(
            passive_tracer   = ds.dim["NX"] - 2,
            tracer_forcings  = tracer_forcings,
            Ny               = ds.dim["Ny"],
            Nz               = ds.dim["Nz"],
            Ω                = ds.attrib["Omega"],
            ϵ                = ds.attrib["epsilon"],
            Δλb              = ds.attrib["dlambda_b"],
            Δλi              = ds.attrib["dlambda_i"],
            ϕn               = ds.attrib["phi_n"],
            ϕs               = ds.attrib["phi_s"],
            H                = ds.attrib["H"],
            R                = ds.attrib["R"],
            Δt               = ds.attrib["dt"],
            Kh               = ds.attrib["Kh"],
            Kv               = ds.attrib["Kv"],
            Kv_cva           = ds.attrib["Kv_cva"],
            cva_Δ            = ds.attrib["cva_delta"],
            ϕc               = ds.attrib["phi_c"],
            Δϕ_trans         = ds.attrib["dphi_trans"],
            Q                = ds.attrib["Q"],
            Q_shape          = ds.attrib["Q_shape"],
            ξ                = ds.attrib["xi"],
            MLT_T            = ds.attrib["MLT_T"],
            MLT_S            = ds.attrib["MLT_S"],
            MLT_shape        = ds.attrib["MLT_shape"],
            μ                = ds.attrib["mu"],
        )

        m = Model(env)

        @fast_extract m

        st.X_ .= ds["X_"]
        st.X_tgt_ .= ds["X_tgt_"]

        updateB!(m)
        updateSaltTransForcing!(m)
    end

    return m
end
