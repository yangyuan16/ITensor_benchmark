using ITensors
using ITensors.HDF5
#
function charge(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    workpath = "./output/psi/"
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_txy$(txy)_tz$(tz)_Jxy$(Jxy)_Jz$(Jz)_dop$(dop_level)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    #
    f = h5open(filename_psi, "r")
    psi = read(f, "psi", MPS)
    close(f)
    #--- measure local electron density 
    Nup = expect(psi, "Nup")
    @show Nup
    Ndn = expect(psi, "Ndn")
    @show Ndn
    Ntot = expect(psi, "Ntot")
    @show Ntot
    #
    savepath = "./"
    filename_Nup_1 = filename_psi_1
    filename_Nup_2 = filename_psi_2
    filename_Nup_3 = "_Nup.h5"
    filename_Nup = join([savepath, filename_Nup_1, filename_Nup_2, filename_Nup_3])

    f = h5open(filename_Nup, "w")
    write(f, "Nup", Nup)
    close(f)
    #----------------------------------------------
    savepath = "./"
    filename_Ndn_1 = filename_psi_1
    filename_Ndn_2 = filename_psi_2
    filename_Ndn_3 = "_Ndn.h5"
    filename_Ndn = join([savepath, filename_Ndn_1, filename_Ndn_2, filename_Ndn_3])

    f = h5open(filename_Ndn, "w")
    write(f, "Ndn", Ndn)
    close(f)
    #----------------------------------------------
    savepath = "./"
    filename_Ntot_1 = filename_psi_1
    filename_Ntot_2 = filename_psi_2
    filename_Ntot_3 = "_Ntot.h5"
    filename_Ntot = join([savepath, filename_Ntot_1, filename_Ntot_2, filename_Ntot_3])

    f = h5open(filename_Ntot, "w")
    write(f, "Ntot", Ntot)
    close(f)
    #
end

#=
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
    println("---psi_obs3_charge---")
    charge(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
end
=#