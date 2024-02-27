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
    #--- measure Sz correlation ---
    corr_z = correlation_matrix(psi, "Sz", "Sz")
    corr_spsm = correlation_matrix(psi, "S+", "S-")
    corr_smsp = correlation_matrix(psi, "S-", "S+")
    #@show corr_z
    #@show corr_spsm
    #@show corr_smsp
    
    # save corrz
    savepath = "./"
    filename_corrz_1 = filename_psi_1
    filename_corrz_2 = filename_psi_2
    filename_corrz_3 = "_corrz.h5"
    filename_corrz = join([savepath, filename_corrz_1, filename_corrz_2, filename_corrz_3])

    f = h5open(filename_corrz, "w")
    write(f, "corrz", corr_z)
    close(f)

    # save corrspsm
    savepath = "./"
    filename_corrspsm_1 = filename_psi_1
    filename_corrspsm_2 = filename_psi_2
    filename_corrspsm_3 = "_corrspsm.h5"
    filename_corrspsm = join([savepath, filename_corrspsm_1, filename_corrspsm_2, filename_corrspsm_3])
    
    f = h5open(filename_corrspsm, "w")
    write(f, "corrspsm", corr_spsm)
    close(f)

    #save corrsmsp
    savepath = "./"
    filename_corrsmsp_1 = filename_psi_1
    filename_corrsmsp_2 = filename_psi_2
    filename_corrsmsp_3 = "_corrsmsp.h5"
    filename_corrsmsp = join([savepath, filename_corrsmsp_1, filename_corrsmsp_2, filename_corrsmsp_3])
    
    f = h5open(filename_corrsmsp, "w")
    write(f, "corrsmsp", corr_smsp)
    close(f)

end