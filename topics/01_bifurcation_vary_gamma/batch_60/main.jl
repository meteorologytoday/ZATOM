using Formatting
using Distributed

label = "wid60"
shapes = [ "balanced_tanh", ]
scan_dirs = [ "pos", ]
ncpu=1
nthread=66

#ξs = [ -10, -9, -8 ,-7 ,-6, -5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.25, -1.0, -0.5, 0.0, 0.5, ]
ξs = [ -5, -4.8, -4.6, -4.4, -4.2, -4.0, 0.0]
#append!(ξs, - collect(0.55:0.05:1.5))
#push!(ξs, 2.0)

ds_exp_lb = 8

@sync for shape in shapes
    for scan_dir in scan_dirs
    for ξ in ξs
    
        println("Doing case ξ = $ξ , scan_dir = $scan_dir")
        #@async run(`julia generate_submit_file.jl
        run(`julia generate_submit_file.jl
            --label     $label
            --p-name    xi
            --p-val     $ξ
            --p-fac     100
            --ds-exp-lb $ds_exp_lb
            --scan-dir  $scan_dir
            --shape     $shape
            --ncpu      $ncpu
            --nthread   $nthread
        `)

    end
    end
end
