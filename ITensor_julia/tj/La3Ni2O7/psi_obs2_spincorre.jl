using ITensors
using ITensors.HDF5
#
function get_bonds(ref1,r,Lz, Ly)
    b = fill(0,(r,2))
    for i = 1:r
        b[i,1] = ref1
        b[i,2] = ref1 + i * Lz * Ly
    end
    println("corre. bonds are:")
    @show b
    return b
end
#
function spincorre(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    workpath = "./output/psi/"
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_txy$(txy)_tz$(tz)_Jxy$(Jxy)_Jz$(Jz)_dop$(dop_level)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    #
    f = h5open(filename_psi, "r")
    psi = read(f, "psi", MPS)
    close(f)
    #--- measure all Sz correlation ---
    #corr_z = correlation_matrix(psi, "Sz", "Sz") # measure all two site correlations
    #corr_spsm = correlation_matrix(psi, "S+", "S-") # measure all two site correlations
    #corr_smsp = correlation_matrix(psi, "S-", "S+") # measure all two site correlations
    #@show corr_z
    #@show corr_spsm
    #@show corr_smsp
    #--- measure Sz correlation along a specific path
    ref_site = 2 # change this parameter !!!
    r = 6 # chang this parameter !!! correlation length
    bonds = get_bonds(ref_site,r,Lz,Ly)        
    corr_z = fill(0.0,(r,3))
    corr_spsm = fill(0.0,(r,3))
    corr_smsp = fill(0.0,(r,3))
    sites = siteinds(psi)
    for j = 1:size(bonds)[1]
        SzSz = AutoMPO()
        SpSm = AutoMPO()
        SmSp = AutoMPO()
        SzSz += 1.0, "Sz", bonds[j,1], "Sz", bonds[j,2]
        SpSm += 1.0, "S+", bonds[j,1], "S-", bonds[j,2]
        SmSp += 1.0, "S-", bonds[j,1], "S+", bonds[j,2]
        #
        SzSz = MPO(SzSz, sites)
        SpSm = MPO(SpSm, sites)
        SmSp = MPO(SmSp, sites)
        #
        corr_z[j,1] = bonds[j,1]
        corr_z[j,2] = bonds[j,2]
        corr_z[j,3] = inner(psi',SzSz,psi)
        #
        corr_spsm[j,1] = bonds[j,1]
        corr_spsm[j,2] = bonds[j,2]
        corr_spsm[j,3] = inner(psi',SpSm,psi)
        #
        corr_smsp[j,1] = bonds[j,1]
        corr_smsp[j,2] = bonds[j,2]
        corr_smsp[j,3] = inner(psi',SmSp,psi)

    end
    @show corr_z
    @show corr_spsm
    @show corr_smsp
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
    println("---psi_obs2_spincorre---")
    spincorre(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
end
=#