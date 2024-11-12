mutable struct Env

    # Number of tracers. NX ≥ 2. First is temperature, second is salinity.
    # Starting the third is the passive tracer
    NX :: Integer 


    Ny_bsn :: Int64
    Nz     :: Int64
    
    Δλb   :: Float64
    Δλi   :: Float64

    gd_b    :: Grid
    gd_bb   :: Grid
    gd_bib  :: Grid

    Δt  :: Float64
    
    Kh      :: Float64
    Kv      :: Float64
    Kv_cva  :: Float64
    Khx     :: AbstractArray{Float64, 2}

    cva_Δ   :: Float64

    ϕc          :: Float64
    ϕc_idx      :: Integer
    Δϕ_trans    :: Float64
    Q           :: Float64
    Q_shape     :: String
    ξ           :: Float64

    MLT_T       :: Float64
    MLT_S       :: Float64
    MLT_shape   :: String

    tracer_forcings :: Array{TracerForcing}

    μ :: Float64 # for EllipticSolver

    function Env(;
        passive_tracer   :: Integer = 0,
        tracer_forcings  :: Array{TracerForcing},
        Ny    :: Int64,
        Nz    :: Int64,
        Ω     :: Float64,
        ϵ     :: Float64,
        Δλb   :: Float64,
        Δλi   :: Float64,
        ϕn    :: Float64, 
        ϕs    :: Float64,
        H     :: Float64,
        R     :: Float64,

        Δt     :: Float64,
        Kh     :: Float64,
        Kv     :: Float64,
        Kv_cva :: Float64,
        Khx :: Union{Nothing, AbstractArray{Float64}} = nothing,
        cva_Δ  :: Float64,
        ϕc       :: Float64,
        Δϕ_trans :: Float64,
        Q        :: Float64,
        Q_shape  :: String,
        ξ        :: Float64,
        MLT_T    :: Float64,
        MLT_S    :: Float64,
        MLT_shape:: String,
        μ        :: Float64 = 1.0,
    )

        s = 0.6
   
        if Δλb <= 0
            throw(ErrorException("`Δλb` must be > 0"))
        end

        if Δλi < 0
            throw(ErrorException("`Δλi` must be >= 0. We got $(Δλi)"))
        end
     
        gd_b = ZAOM.Grid(;
            ϵ   = ϵ,
            Ω   = Ω,
            Ny  = Ny,
            Nz  = Nz,
            Δλ  = [Δλb,],
            ϕn  = ϕn,
            ϕs  = ϕs,
            H   = H,
            R   = R,
            s   = s,
        )

        gd_bb = ZAOM.Grid(;
            ϵ   = ϵ,
            Ω   = Ω,
            Ny  = Ny,
            Nz  = Nz,
            Δλ  = [Δλb, Δλb],
            ϕn  = ϕn,
            ϕs  = ϕs,
            H   = H,
            R   = R,
            s   = s,
        )

        gd_bib = ZAOM.Grid(;
            ϵ   = ϵ,
            Ω   = Ω,
            Ny  = Ny,
            Nz  = Nz,
            Δλ  = [Δλb, Δλb+Δλi],
            ϕn  = ϕn,
            ϕs  = ϕs,
            H   = H,
            R   = R,
            s   = s,
        )





        # Zonal diffusivity
        if Khx == nothing 

            Khx = zeros(Float64, Nz, Ny)

        elseif length(size(Khx)) == 1

            if length(Khx) != Ny
                throw(ErrorException("Wrong Khx dimension."))
            end

            Khx = repeat( reshape(Khx, 1, Ny), dims=(Nz, 1) )

        elseif length(size(Khx)) == 2

            if size(Khx) != (Nz, Ny)
                throw(ErrorException("Wrong Khx dimension."))
            end

        else
            throw(ErrorException("Unknown Khx specification."))
        end
 
        # North and south boundary to force satisfying
        # the boundary condition be == bw
        #Khx[:,   :]   .= 5e4 * exp.( - ( gd_b.ϕ_T[:, :, 1] / deg2rad(20) ).^2 / 2 )
        Khx[:,   1]   .= 1e10
        Khx[:, end]   .= 1e10

        NX = passive_tracer + 2

        if length(tracer_forcings) != NX
            throw(ErrorException("TracerForcing array should have " * string(NX) * " members."))
        end

        # Round-up ϕc
        ϕ_T = gd_bib.ϕ_T[1, :, 1]
        ϕ_V = gd_bib.ϕ_V[1, :, 1]
        ϕc_idx = findfirst(ϕc .< ϕ_T) - 1
        println(format("ϕc_idx = {:d}. So the actual ϕc changes from {:.1f} to {:.1f} degrees", ϕc_idx, rad2deg(ϕc), rad2deg(ϕ_V[ϕc_idx+1])))
        ϕc = ϕ_V[ϕc_idx + 1]

        return new(
            NX,
            Ny,
            Nz,
            Δλb,
            Δλi,
            gd_b,
            gd_bb,
            gd_bib,
            Δt,
            Kh,
            Kv,
            Kv_cva,
            Khx,
            cva_Δ,
            ϕc,
            ϕc_idx,
            Δϕ_trans,
            Q,
            Q_shape,
            ξ,
            MLT_T,
            MLT_S,
            MLT_shape,
            tracer_forcings,
            μ,
        )   
    end
end

