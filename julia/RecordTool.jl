module RecordTool

    using NCDatasets
    using Formatting
    using Dates
    missing_value = 1e20
 
    mutable struct VarInfo
        varname  :: AbstractString
        dimnames :: Union{Nothing, Tuple}   # nothing = scalar
        attr     :: Union{Nothing, Dict}   
    end

    mutable struct VarObj

        vartype  :: Symbol    # :STAT :NONSTAT
        varinfo  :: VarInfo
        varref   :: AbstractArray{Float64}
        var      :: AbstractArray{Float64}
        weight   :: Float64

        function VarObj(
            vartype :: Symbol,
            varinfo :: VarInfo,
            varref  :: AbstractArray{Float64},
        )
            if ! ( vartype  in ( :STAT, :NONSTAT ) )
                throw(ErrorException("vartype only allows :STAT and :NONSTAT"))
            end

            if length(size(varref)) != length(varinfo.dimnames)
                ErrorException(format("Variable `{:s}` is {:d} dimensional while only {:d} dimension names are given.", varinfo.varname, length(size(varref)), length(varinfo.dimnames))) |> throw
            end

            var = zeros(Float64, size(varref)...)

            return new(vartype, varinfo, varref, var, 0.0)
        end
    
    end

    
    mutable struct Recorder
        
        dict_varinfos :: Dict

        filename :: Union{Nothing, AbstractString}
        time_ptr :: Integer  # The position of next record
        
        dims     :: Dict    # A dictionary of dimension name mapping to its length
        sobjs    :: Dict
        nsobjs   :: Dict

        function Recorder(dims, varinfos; other_vars = nothing)

            sobjs  = Dict()
            nsobjs = Dict()

            if haskey(dims, "time")
                ErrorException("Dimension `time` is used for record dimension. It cannot be specified in dict `dims`.") |> throw
            end

            dict_varinfos = Dict()
            for varinfo in varinfos
                dict_varinfos[varinfo.varname] = varinfo
            end

            return new(dict_varinfos, nothing, 1, dims, sobjs, nsobjs)

        end



    end


    function bindVariable!(
        rec      :: Recorder,
        vartype  :: Symbol,
        varname  :: AbstractString,
        varref   :: AbstractArray{Float64},
    )

        for dimname in varinfo.dimnames
            if !haskey(dims, dimname)
                ErrorException(format("Variable `{:s}` contains a dimension `{:s}` not specified in `dims` dict.", varinfo.varname, dimname)) |> throw
            end
        end

        rec.sobjs[varname] = VarObj(vartype, rec.dict_varinfos[varname], varref)
    end

    function bindVariables!(
        rec      :: Recorder;
        stats    :: Union{Array{Dict}, Nothing} = nothing,
        nonstats :: Union{Array{Dict}, Nothing} = nothing,
    )

        if stats != nothing
            for stat in stats
                bindVariable!(rec, :STAT, stat.first, stat.second)
            end
        end

        if nonstats != nothing
            for nonstat in nonstats
                bindVariable!(rec, :NONSTAT, nonstat.first, nonstat.second)
            end
        end
    
    end

    function setNewNCFile!(
        rec      :: Recorder,
        filename :: AbstractString,
    )
        
        rec.filename = filename
        rec.time_ptr = 1
        
        Dataset(filename, "c") do ds
        
            for (dimname, dim) in rec.dims 
                defDim(ds, dimname, dim)
            end

            defDim(ds, "time", Inf)
            ds.attrib["_FillValue"] = missing_value
            ds.attrib["timestamp"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS sss")

            for (varname, sobj) in rec.sobjs
                ds_var = defVar(ds, varname, Float64, (sobj.varinfo.dimnames..., "time"))
                ds_var.attrib["_FillValue"] = missing_value

                if sobj.varinfo.attr != nothing
                    for (k, v) in sobj.varinfo.attr
                        ds_var.attrib[k] = v
                    end
                end
            end

            for (varname, nsobj) in rec.nsobjs
                ds_var = defVar(ds, varname, Float64, (nsobj.dimnames...,))
                ds_var.attrib["_FillValue"] = missing_value

                if sobj.varinfo.attr != nothing
                    for (k, v) in sobj.varinfo.attr
                        ds_var.attrib[k] = v
                    end
                end
                
                ds_var[:] = nsobj.varref

            end 

 
        end
        
    end

    function record!(
        rec::Recorder;
    )

        varnames = keys(rec.sobjs)

        for varname in varnames
            sobj = rec.sobjs[varname]
            sobj.var .+= sobj.varref
            sobj.weight += 1.0
        end
        
    end

#=
    function record_wrap!(
        rec             :: Recorder;
        create_new_file :: Bool,
        avg_and_output  :: Bool,
        new_file_name   :: AbstractString,
    )

        if create_new_file
            setNewNCFile!(rec, new_file_name)
        end

        record!(rec; avg_and_output=avg_and_output)

    end
=#
 

    function avgAndOutput!(
        rec :: Recorder
    )

        if rec.filename == nothing
            ErrorException("Undefined record filename") |> throw
        end

        # Do average
        for (varname, sobj) in rec.sobjs
            if sobj.weight == 0
                ErrorException(format("StatObj for variable `{:s}` has weight 0 during normalization.", varname)) |> throw
            end
            sobj.var /= sobj.weight
        end
        
        # Output data
        Dataset(rec.filename, "a") do ds
            for (varname, sobj) in rec.sobjs
                ds_var = ds[varname]
                ds[varname][repeat([:,], length(sobj.dimnames))..., rec.time_ptr] = sobj.var
            end
        end


        # Reset StatObjs
        for (_, sobj) in rec.sobjs
            sobj.var .= 0.0
            sobj.weight = 0.0
        end
        
        # Increment of time
        rec.time_ptr += 1
    
    end

end
