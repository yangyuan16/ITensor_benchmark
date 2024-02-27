using ITensors
using ITensors.HDF5
#
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
    workpath = "./output/psi/"
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_txy$(txy)_tz$(tz)_Jxy$(Jxy)_Jz$(Jz)_dop$(dop_level)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    #
    f = h5open(filename_psi, "r")
    psi = read(f, "psi", MPS)
    close(f)
    #--- measure charge density correlation function
    charge_density_corre = correlation_matrix(psi, "Ntot", "Ntot")
    #@show Green_up
    #@show Green_dn
    
    savepath = "./"
    filename_chargedensitycorre_1 = filename_psi_1
    filename_chargedensitycorre_2 = filename_psi_2
    filename_chargedensitycorre_3 = "_Chargedensitycorre.h5"
    filename_chargedensitycorre = join([savepath, filename_chargedensitycorre_1, 
                              filename_chargedensitycorre_2, filename_chargedensitycorre_3])

    f = h5open(filename_chargedensitycorre, "w")
    write(f, "chargedensitycorre", charge_density_corre)
    close(f)
end