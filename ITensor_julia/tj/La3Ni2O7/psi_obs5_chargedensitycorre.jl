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
function chargedensitycorre(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
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
    #charge_density_corre = correlation_matrix(psi, "Ntot", "Ntot") # measure all two sites density correlation  
    #@show Green_up
    #@show Green_dn
    #--- measure Green function along a specific path
    ref_site = 2 # change this parameter !!!
    r = 6 # chang this parameter !!! correlation length
    bonds = get_bonds(ref_site,r,Lz,Ly)        
    charge_density_corre = fill(0.0,(r,3))
    sites = siteinds(psi)
    for j = 1:size(bonds)[1]
        D_corr = AutoMPO()
        D_corr += 1.0, "Ntot", bonds[j,1], "Ntot", bonds[j,2]
        #
        D_corr = MPO(D_corr, sites)
        #
        charge_density_corre[j,1] = bonds[j,1]
        charge_density_corre[j,2] = bonds[j,2]
        charge_density_corre[j,3] = inner(psi',D_corr,psi)
    end
    @show charge_density_corre
    #
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
    println("---psi_obs5_chargedensitycorre---")
    chargedensitycorre(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
end
=#