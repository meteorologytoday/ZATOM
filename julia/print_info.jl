function printInfo(m::Model)
    @fast_extract m

    println("===== Model Info =====")
    println(format("Domain (Nz, Ny) = ({:d}, {:d})", ev.gd.Nz, ev.gd.Ny))

end
