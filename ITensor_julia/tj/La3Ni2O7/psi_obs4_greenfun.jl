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
function greenfun(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
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
        #Green_up = correlation_matrix(psi, "Cdagup", "Cup") # measure all two site single particle green function
        #Green_dn = correlation_matrix(psi, "Cdagdn", "Cdn") # measure all two site single particle green function
        #@show Green_up
        #@show Green_dn
        #--- measure Green function along a specific path
        ref_site = 2 # change this parameter !!!
        r = 6 # chang this parameter !!! correlation length
        bonds = get_bonds(ref_site,r,Lz,Ly)        
        Green_up = fill(0.0,(r,3))
        Green_dn = fill(0.0,(r,3))
        sites = siteinds(psi)
        for j = 1:size(bonds)[1]
            Gup = AutoMPO()
            Gdn = AutoMPO()
            Gup += 1.0, "Cdagup", bonds[j,1], "Cup", bonds[j,2]
            Gdn += 1.0, "Cdagdn", bonds[j,1], "Cdn", bonds[j,2]
            #
            Gup = MPO(Gup, sites)
            Gdn = MPO(Gdn, sites)
            #
            Green_up[j,1] = bonds[j,1]
            Green_up[j,2] = bonds[j,2]
            Green_up[j,3] = inner(psi',Gup,psi)
            #
            Green_dn[j,1] = bonds[j,1]
            Green_dn[j,2] = bonds[j,2]
            Green_dn[j,3] = inner(psi',Gdn,psi)
            #
        end
        @show Green_up
        @show Green_dn
        #
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
    println("---psi_obs4_greenfun---")
    greenfun(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
end
=#