function genSaltTransForcing(
    env :: Env;
)
    tol = 1e-5 # tolerance for freshwater flux balance

    gd = env.gd_bib
    MLT_w = genMixedlayerWeight(gd.z_T, env.MLT_S, env.MLT_shape) 


    ϕ_T = gd.ϕ_T

    # Mask is assumed to be the whole ocean
    mask = ones(Float64, size(ϕ_T)...)

    # Surface mask includes the grid points within H_sfc
    surf_mask = zeros(Float64, size(ϕ_T)...)
    surf_mask[MLT_w .!= 0] .= 1.0

    ϕn = maximum(gd.ϕ_T)
    ϕs = minimum(gd.ϕ_T)
        
    σ = similar(ϕ_T)
    σ .= 0.0

    balanced_mode, γ_shape = split(env.Q_shape, "_")
    
    if γ_shape == "flat"
        σ[:, (env.ϕc_idx+1):end, 1] .= 1.0
    elseif γ_shape == "linear"
        σ .= ϕ_T
    elseif γ_shape == "cosine"
        println(rad2deg.(ϕ_T[1, :, 1]))
        σ .= - cos.(π * (ϕ_T .- ϕs) / (ϕn - ϕs) )
    elseif γ_shape == "tanh"
        σ .= tanh.((ϕ_T .- env.ϕc) / env.Δϕ_trans)
    else
        throw(ErrorException("Unrecognized γ_shape: $γ_shape."))
    end

    # remove the σ below the mixed layer
    @. σ *= MLT_w * surf_mask

    if balanced_mode == "raw"
        # do nothing
    elseif balanced_mode == "balanced"

        # Remember σ is already filtered by surf_mask so that
        # we do not have to worry about the grids below H_sfc
 
        σΔv_T = σ .* gd.Δv_T
        ∫_posσΔv_T = sum(σΔv_T[σ .>= 0])
        ∫_negσΔv_T = sum(σΔv_T[σ .< 0])

        if ∫_posσΔv_T > 0
            #println("Doing adjustment")
            #println("∫_posσΔv_T = $∫_posσΔv_T")
            #println("∫_negσΔv_T = $∫_negσΔv_T")
            factor = - ∫_negσΔv_T / ∫_posσΔv_T
            println("Balance factor: $factor")
            σ[σ .>= 0 ] .*= factor
        end

        # Verify
        σΔv_T = σ .* gd.Δv_T
        ∫_posσΔv_T = sum(σΔv_T[σ .>= 0])
        ∫_negσΔv_T = sum(σΔv_T[σ .< 0])   # Redundant but I keep it for clarity
        ∫_σΔv_T_surf = ∫_posσΔv_T + ∫_negσΔv_T
        err =  abs( ( ∫_σΔv_T_surf - 0 ) / ∫_posσΔv_T)  # Measure in terms of ratio
        if ∫_posσΔv_T > 0 && err > tol
            println("∫_posσΔv_T = $∫_posσΔv_T")
            println("∫_negσΔv_T = $∫_negσΔv_T")
            println("Error: abs( (∫_σΔv_T_surf - 0) / ∫_posσΔv_T) = $err > $tol")
            throw(ErrorException("Balanced mode is on, but I cannot balance the forcing."))
        end
    else
        throw(ErrorException("Unrecognized balanced mode: $balanced_mode."))
    end

    # ===== Adjust so that ∫ σ dv = 0 =====
    σ̅ = sum(mask .* gd.Δv_T .* σ) / sum(mask .* gd.Δv_T)
    σ[mask .== 1] .-= σ̅

    # Computing r factor 
    # Adjust east-west ratio
    r = similar(σ)
    drdξ = similar(σ)

    # Remember to apply surf_mask. We do not want to
    # alter the flux ratio in deep ocean.
    west_pos = (view(σ, :, :, 1) .>= 0.0) .& ( view(surf_mask, :, :, 1) .== 1.0)
    east_pos = (view(σ, :, :, 2) .>= 0.0) .& ( view(surf_mask, :, :, 2) .== 1.0)

    D = env.Δλb / (env.Δλi + env.Δλb)
    r .= 1.0
    view(r, :, :, 1)[west_pos] .= 1.0 - env.ξ
    view(r, :, :, 2)[east_pos] .= 1.0 + env.ξ * D

    drdξ .= 0.0
    view(drdξ, :, :, 1)[west_pos] .= - 1
    view(drdξ, :, :, 2)[east_pos] .=   D

    
    # ===== Adjust so that ∫ σ dv (integrate over σ > 0) = 1.0 =====
    σΔv_T = σ .* r .* gd.Δv_T
    ∫_posσΔv_T = sum(σΔv_T[σ .>= 0])
    if ∫_posσΔv_T != 0
        σ ./= ∫_posσΔv_T
    end
 
    # Computing resulting matrix
    Sf    = zeros(Float64, gd.Nz, gd.Ny, gd.Nx)
    dSfdQ = copy(Sf)
    dSfdξ = copy(Sf)
   
    # For clarity, let's compute it in a tedious way
    @. Sf    = - S0 * σ * r  * env.Q
    @. dSfdQ = - S0 * σ * r
    @. dSfdξ = - S0 * σ * drdξ * env.Q
 
    # ===== Verify ===== 
    qΔv_T = env.Q * σ .* r .* gd.Δv_T
    ∫_posqΔv_T = sum(qΔv_T[qΔv_T .>= 0])
    ∫_qΔv_T = sum(qΔv_T)

    err1  = abs( ( ∫_posqΔv_T - env.Q ) / env.Q ) 
    err2 = abs( (∫_qΔv_T - 0 ) / ∫_posqΔv_T )
    
    if err1 > tol || ( ∫_posqΔv_T > 0 && err2 > tol )
        println("Warning: The constructed total freshwater flux does not match with target.")
        println("env.Q = ", env.Q) 
        println("∫_posqΔv_T = ", ∫_posqΔv_T) 
        println("∫_qΔv_T = ", ∫_qΔv_T) 
        println("abs(( ∫_posqΔv_T - 0) / env.Q     ) = ", err1) 
        println("abs(( ∫_qΔv_T    - 0) / ∫_posqΔv_T) = ", err2) 
        throw(ErrorException("Error : the specified Q not achieved."))
    end
    # ==================
 
    return Sf, dSfdQ, dSfdξ

end


function updateSaltTransForcing!(
    m  :: Model;
    Q  :: Union{Float64, Nothing} = nothing,
    ξ  :: Union{Float64, Nothing} = nothing, 
)

    if Q != nothing
        m.env.Q = Q
    end

    if ξ != nothing
        m.env.ξ = ξ
    end


    Sf, dSfdQ, dSfdξ = genSaltTransForcing(m.env)
    m.state.X_SS_[:, :, :, 2] .= Sf
    m.core.op_mtx[:dSfdQ] = dSfdQ
    m.core.op_mtx[:dSfdξ] = dSfdξ
 
end


function calTotalFreshwaterFlux(
    m :: Model,
)
    salt_forcing = view(m.state.X_SS_, :, :, :, 2)
    fwf = - salt_forcing / S0

    idx = fwf .> 0.0
    fwf = sum(m.env.gd_bib.Δv_T[idx] .* fwf[idx])
    return fwf
end


