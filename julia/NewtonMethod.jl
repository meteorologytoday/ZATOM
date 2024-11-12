module NewtonMethod

    export doNewtonMethodOnce!, doNewtonMethod!

    using LinearAlgebra: logabsdet, norm
    using Formatting


    function doNewtonMethodOnce!(
        X    :: Array,
        F_dFdX :: Function;
        callback :: Union{Function, Nothing} = nothing,
        Δ    :: Any = nothing,
        J    :: Any = nothing,
        verbose :: Bool = false,
    )


        # This design is to avoid re-compute Δ when
        # doing multiple newton steps in a row. The
        # Δ used next time is the Δ after stepped 
        # this time.
        if Δ == nothing || J == nothing
            Δ, J = F_dFdX(X)
        end

        X .-= J \ Δ

        if callback != nothing
            callback()
        end


        Δ, J = F_dFdX(X)

        res = norm(Δ) / length(Δ)

        detJ_info = logabsdet(Matrix(J))

        if verbose
            println(format("res={:e}, log(D)={:e}", res, detJ_info[1]))
        end


        return res, Δ, J, detJ_info

    end

    function doNewtonMethod!(
        X          :: Array,
        F_dFdX     :: Function,
        res_target :: Float64,
        min_iter   :: Integer,
        max_iter   :: Integer;
        callback   :: Union{Function, Nothing} = nothing,
        verbose        :: Bool = false,
    )

        if max_iter < min_iter
            throw(ErrorException("`max_iter` must be greater or equal to `min_iter`."))
        end

        cnt = 0
        res = Inf
        detJ_info = nothing

        flag = nothing

        Δ = nothing
        J = nothing

        while true
     
            if cnt >= min_iter && res <= res_target
                
                verbose && println(format("Converge, now breaking... cnt={:d}, res={:e}", cnt, res))
                flag = :converge
                break
                 
            end
           
            if cnt < max_iter

                res, Δ, J, detJ_info = doNewtonMethodOnce!(X, F_dFdX; Δ=Δ, J=J, verbose=verbose, callback=callback)
                if callback != nothing
                    callback()
                end
                cnt += 1

            else

                flag = :no_converge
                break

            end

            if !isfinite(res) 
                
                flag = :no_converge
                break

            end

        end
        
        return flag, cnt, res, Δ, J, detJ_info

    end

end
