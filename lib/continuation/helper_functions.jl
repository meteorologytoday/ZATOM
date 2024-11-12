
function detectInitialization(
    wdir     :: String,
    casename :: String,
)

    local max_cnt  = -1
    local filename = nothing
    local next_filename = nothing
    local status   = nothing
    local casedir = joinpath(wdir, casename)
    if !isdir(casedir)
        mkdir(casedir)
    end
    for _filename in readdir(casedir)
        println(_filename) 
        local m = match(r"^(?<num>[0-9]+)\.nc$", _filename)
        if m != nothing
            num = parse(Int64, m[:num])
            println("Find $(_filename) whose count is $(num).")
            
            if max_cnt < num
                max_cnt = num
                filename = joinpath(casedir, _filename)
                next_filename = joinpath(casedir, format("{:04d}.nc", max_cnt+1))
            end
        end

    end

    if max_cnt == -1
        status = :INIT
        next_filename = joinpath(casedir, format("{:04d}.nc", 0))
    elseif max_cnt == 9999
        throw(ErrorException("Max cnt is 9999 which is the maximum."))
    else
        status = :CONT
    end

    return status, filename, next_filename

end


