using Formatting
using Distributed

label = "Xi_60"
shapes = [ "balanced_tanh", ]
scan_dirs = [ "neg", ]
ncpu=1
nthread=66


γs = collect(0:0.025:0.2)[2:end]
push!(γs, 0.005)


ds_exp_lb = 8

@sync for shape in shapes
    for scan_dir in scan_dirs
    for γ in γs
    
        println("Doing case γ = $γ , scan_dir = $scan_dir")
        #@async run(`julia generate_submit_file.jl
        run(`julia generate_submit_file.jl
            --gamma     $γ
            --ds-exp-lb $ds_exp_lb
            --scan-dir  $scan_dir
            --shape     $shape
            --ncpu      $ncpu
            --nthread   $nthread
        `)

    end
    end
end


# positive cases
scan_dirs = [ "pos", ]

γs = [0.005, 0.025]


ds_exp_lb = 8

@sync for shape in shapes
    for scan_dir in scan_dirs
    for γ in γs
    
        println("Doing case γ = $γ , scan_dir = $scan_dir")
        #@async run(`julia generate_submit_file.jl
        run(`julia generate_submit_file.jl
            --gamma     $γ
            --ds-exp-lb $ds_exp_lb
            --scan-dir  $scan_dir
            --shape     $shape
            --ncpu      $ncpu
            --nthread   $nthread
        `)

    end
    end
end
