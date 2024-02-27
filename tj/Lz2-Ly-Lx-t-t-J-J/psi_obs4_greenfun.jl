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
    #--- measure single particle Green function
    Green_up = correlation_matrix(psi, "Cdagup", "Cup")
    Green_dn = correlation_matrix(psi, "Cdagdn", "Cdn")
    #@show Green_up
    #@show Green_dn
    
    savepath = "./"
    filename_Greenup_1 = filename_psi_1
    filename_Greenup_2 = filename_psi_2
    filename_Greenup_3 = "_Greenup.h5"
    filename_Greenup = join([savepath, filename_Greenup_1, filename_Greenup_2, filename_Greenup_3])

    f = h5open(filename_Greenup, "w")
    write(f, "Green_up", Green_up)
    close(f)

    savepath = "./"
    filename_Greendn_1 = filename_psi_1
    filename_Greendn_2 = filename_psi_2
    filename_Greendn_3 = "_Greendn.h5"
    filename_Greendn = join([savepath, filename_Greendn_1, filename_Greendn_2, filename_Greendn_3])

    f = h5open(filename_Greendn, "w")
    write(f, "Green_dn", Green_dn)
    close(f)


end