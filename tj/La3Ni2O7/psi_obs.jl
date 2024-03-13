using ITensors
using ITensors.HDF5
#
include("psi_info.jl")
include("psi_obs1_magz.jl")
include("psi_obs2_spincorre.jl")
include("psi_obs3_charge.jl")
include("psi_obs4_greenfun.jl")
include("psi_obs5_chargedensitycorre.jl")
include("psi_obs6_singletpairing.jl")
let 
    # load psi
    bc = "pbc"
    Lz = 2
    Ly = 3
    Lx = 8
    N = Lz * Ly * Lx
    #
    txy = 3.0
    tz = 0.0
    Jxy = 1.0
    Jz = 0.5
    dop_level = (1,2)
    #
    println("---psi_info---")
    info(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    println("---psi_obs1_magz---")
    magz(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    println("---psi_obs2_spincorre---")
    spincorre(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    println("---psi_obs3_charge---")
    charge(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    println("---psi_obs4_greenfun---")
    greenfun(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    println("---psi_obs5_chargedensitycorre---")
    chargedensitycorre(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    println("---psi_obs6_singletpairing---")
    singletpairing(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    
end