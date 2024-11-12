Nz = 20
Ny = 30

cfg = SortedDict{String, Any}(
    "ϵ"           => 0.0,#1e-4,
    "Nz"          => Nz,
    "Ny"          => Ny,
    "MLT_T"       => 60.0,
    "MLT_S"       => 60.0,
    "MLT_shape"   => "step",
    "domain_size" => (basin=50.0, bnd=5.0),
    "γ_rng_pos"   => [0.0, 0.25],
    "γ_rng_neg"   => [0.0, 2.0],
    "lat"         => [10.0, 70.0],
    "ϕc"          => 40.0,
    "Δϕ_trans"    => 10,
    "Kh"          => 4e4,
    "Kv_iso"      => 1e-4,
    "spinup"      => Dict("pos" => 3, "neg" => 250), # years
    "ξ0"          => -1.0,
    "γ0"          =>  0.0,
    "scales"      => [ 1.0, 1.0, 0.01e6 ],
    "cva_Δ"       => 1e-4,
    "use_transient_ξ_pos" => false,
    "use_transient_ξ_neg" => false,
    "transient_ξ" => 2.0,
    "transient_ξ_duration" => 50,
    "T_restore_days" => 5,
    "dT-east"     => 25.0,
    "dT-west"     => 25.0,
    "μ"           => 1.0,
)
