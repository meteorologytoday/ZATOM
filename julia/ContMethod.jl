if :NewtonMethod ∉ names(Main)
    include(joinpath(@__DIR__, "NewtonMethod.jl"))
end


module ContMethod1

    export CMInfo, doContinuition!, setS!, setXP!

    using ..NewtonMethod
    using LinearAlgebra
    using Formatting
    using SparseArrays

    include("WorkspaceContMethod.jl")

    mutable struct CMInfo

        Nx   :: Int64
        
        Mx   :: AbstractArray{Float64, 2}
        Mp   :: AbstractArray{Float64, 2}
       
        invMx :: AbstractArray{Float64, 2}
        invMp :: AbstractArray{Float64, 2}
 
        func     :: Any

        nwt_min_iter :: Int64
        nwt_max_iter :: Int64
        res      :: Float64

        newton_callback :: Union{Function, Nothing}

        s       :: AbstractArray{Float64, 1}
        ṡ       :: AbstractArray{Float64, 1}
        ds_exp  :: Float64    # ds = 2^{ds_exp}
        Δds_dn_exp :: Float64    # The change of ds_exp when decreasing
        Δds_up_exp :: Float64    # The change of ds_exp when increasing
        ds_exp_lb :: Float64  # lower bound of ds_exp

        skip_unconverge :: Bool

        wksp :: Workspace

        function CMInfo(;
            Nx   :: Int64,
            F_dFdx_dFdp :: Function,
            mx   :: AbstractArray{Float64, 1}, 
            mp   :: Float64,
            nwt_min_iter :: Int64,
            nwt_max_iter :: Int64,
            res  :: Float64,
            newton_callback :: Union{Function, Nothing} = nothing,
            Δds_dn_exp   :: Float64 = 1.0,
            Δds_up_exp   :: Float64 = .1,
            ds_exp_lb :: Float64 = -10.0,
            skip_unconverge :: Bool = false,
        )
            if length(mx) != Nx
                throw(ErrorException("Length of mx should equal to Nx."))
            end

            Mx = spdiagm( 0 => mx )
            Mp = spdiagm( 0 => [ mp ] )
 
            invMx = spdiagm( 0 => (1.0 ./ mx) )
            invMp = spdiagm( 0 => [ 1.0 / mp ] )

            s2xp(s) = (Mx*view(s, 1:Nx), Mp*view(s, Nx+1:Nx+1))

            _Λ_dΛds(s, ṡ, s0, ds)  = (
                ṡ' * (s - s0) .- ds,
                ṡ',
            )

            _G_dGds = function(s, ṡ, s0, ds)
           
                local F, dFdx, dFdp = F_dFdx_dFdp(s2xp(s)...)
                local dFds = [ dFdx*Mx  dFdp*Mp ]
                local Λ, dΛds = _Λ_dΛds(s, ṡ, s0, ds) 

                G = [
                    F ;
                    Λ ; 
                ]

                dGds = [
                    dFds ;
                    dΛds ;
                ]

                return G, dGds 
            end

            _ṡ = function(ṡ0, s)
                return normalize(ṡ0)
            end

            func = (
                _Λ_dΛds = _Λ_dΛds,
                _G_dGds = _G_dGds,
                _ṡ    = _ṡ,
            )

            s        = zeros(Float64, Nx + 1)
            ṡ        = zeros(Float64, Nx + 1)
            ds_exp = 0.0
            if Δds_dn_exp <= 0
                throw(ErrorException("Δds_dn_exp should be a positive number."))
            end
            if Δds_up_exp <= 0
                throw(ErrorException("Δds_dn_exp should be a positive number."))
            end
            if ds_exp_lb > 0
                throw(ErrorException("ds_exp_lb should be a negative number."))
            end



            return new(
                Nx,
                Mx,
                Mp,
                invMx,
                invMp,
                func,
                nwt_min_iter,
                nwt_max_iter,
                res,
                newton_callback,

                s,
                ṡ,
                ds_exp,
                Δds_dn_exp,
                Δds_up_exp,
                ds_exp_lb,
                skip_unconverge,
 
                Workspace(Nx+1),
            )

        end
    end

    function setS!(
        cmi :: CMInfo;
        s :: Array{Float64, 1},
        x :: Union{Nothing, Array{Float64, 1}} = nothing,
        p :: Union{Nothing, Array{Float64, 1}} = nothing,
        no_scale :: Bool = false,
    )
        if no_scale
            if x != nothing
                s[1:end-1] .= x
            end
            if p != nothing
                s[end] = p[1]
            end 
        else
            if x != nothing
                mul!(view(s, 1:cmi.Nx), cmi.invMx, x)
            end

            if p != nothing
                mul!(view(s, cmi.Nx+1:cmi.Nx+1), cmi.invMp, p)
            end
        end
    end

    """
        This function sets X and P from s.
    """
    function setXP!(
        cmi :: CMInfo;
        s :: Array{Float64, 1},
        x :: Union{Nothing, Array{Float64, 1}} = nothing,
        p :: Union{Nothing, Array{Float64, 1}} = nothing,
    )
        if x != nothing
            mul!(x, cmi.Mx, view(s, 1:cmi.Nx))
        end

        if p != nothing
            mul!(p, cmi.Mp, view(s, cmi.Nx+1:cmi.Nx+1))
        end
    end



    function doContinuition!(
        cmi :: CMInfo;
        verbose :: Bool = false,
    )

        local dGds
        local converge = false
        local res = -1.0
        local detJ_info
 
        reset!(cmi.wksp)

        # The following two temp fields are necessary
        # because during the callback I want to update
        # cmi.s for some other models might need to 
        # update their internal state to give proper 
        # Jacobian
        s0_tmp = getSpace!(cmi.wksp)
        ṡ0_tmp = getSpace!(cmi.wksp)
   
        normalize!(cmi.ṡ)
        s0_tmp[:] = cmi.s
        ṡ0_tmp[:] = cmi.ṡ

        while true

            local ds = 2.0^cmi.ds_exp 

            if verbose
                println("ds_exp = $(cmi.ds_exp)")
            end

            _G_dGds(s) = cmi.func._G_dGds(s, ṡ0_tmp, s0_tmp, ds)

            # Predictor (Euler forward)
            @. cmi.s += cmi.ṡ * ds

            flag, cnt, res, G_res, dGds, detJ_info = NewtonMethod.doNewtonMethod!(
                cmi.s,
                _G_dGds,
                cmi.res,
                cmi.nwt_min_iter,
                cmi.nwt_max_iter;
                callback = function()
                    cmi.newton_callback(cmi)
                end,
                verbose = verbose,
            )
        
            if flag == :converge
                converge = true
                break
            else

                println("Warning: not converging.")
                if cmi.ds_exp > cmi.ds_exp_lb
                    println("Warning: Change ds_exp from $(cmi.ds_exp) to $(cmi.ds_exp - cmi.Δds_dn_exp).")
                    cmi.ds_exp = max(cmi.ds_exp - cmi.Δds_dn_exp, cmi.ds_exp_lb)

                    # restoring state
                    cmi.s[:] = s0_tmp
                    cmi.ṡ[:] = ṡ0_tmp
                    
                    # remember to callback
                    cmi.newton_callback(cmi)
                else

                    if cmi.skip_unconverge
                        println("Warning: accept this result even it is not converging.")
                        break
                    else
                        throw(ErrorException("Continuition method not converging even ds_exp reaches lower bound."))
                    end

                end

            end
        end

        if converge
            cmi.ds_exp = min(cmi.ds_exp + cmi.Δds_up_exp, 0.0)
        end

        dGds_RHS = getSpace!(cmi.wksp)
        dGds_RHS[1:end-1] .= 0.0
        dGds_RHS[end]      = 1.0

        cmi.ṡ[:] = dGds \ dGds_RHS
        normalize!(cmi.ṡ)

        return converge, res, detJ_info
    end

end
