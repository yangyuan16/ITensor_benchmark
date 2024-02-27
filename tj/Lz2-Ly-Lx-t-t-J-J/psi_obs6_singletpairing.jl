using ITensors
using ITensors.HDF5
#
function get_singletbonds(ref1,ref2,r,Lz, Ly)
    b = fill(0,(r,2))
    for i = 1:r
        b[i,1] = ref1 + i * Lz * Ly
        b[i,2] = ref2 + i * Lz * Ly
    end
    println("singlet bonds are:")
    @show b
    return b
end
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
    #--- measure singlet pairing correlation function
    ref_site1 = 1
    ref_site2 = 2
    r = 6
    singlet_bonds = get_singletbonds(ref_site1,ref_site2,r,Lz,Ly)
    
    pairing = fill(0.0,(r,5))
    sites = siteinds(psi)
    for j = 1:size(singlet_bonds)[1]
        delta = AutoMPO()
        delta += 0.5, "Cdagup", ref_site1, "Cdagdn", ref_site2, "Cup", singlet_bonds[j,1], "Cdn", singlet_bonds[j,2]
        delta += -0.5, "Cdagup", ref_site1, "Cdagdn", ref_site2, "Cdn", singlet_bonds[j,1], "Cup", singlet_bonds[j,2]
        delta += -0.5, "Cdagdn", ref_site1, "Cdagup", ref_site2, "Cup", singlet_bonds[j,1], "Cdn", singlet_bonds[j,2]
        delta += 0.5, "Cdagdn", ref_site1, "Cdagup", ref_site2, "Cdn", singlet_bonds[j,1], "Cup", singlet_bonds[j,2]
        delta = MPO(delta,sites)
        pcorre = inner(psi',delta,psi)
        pairing[j,1] = ref_site1
        pairing[j,2] = ref_site2
        pairing[j,3] = singlet_bonds[j,1]
        pairing[j,4] = singlet_bonds[j,2]
        pairing[j,5] = pcorre 
    end
    @show pairing

    savepath = "./"
    filename_pairing_1 = filename_psi_1
    filename_pairing_2 = filename_psi_2
    filename_pairing_3 = "_pair_ref$((ref_site1,ref_site2)).h5"
    filename_pairing = join([savepath, filename_pairing_1, 
                              filename_pairing_2, filename_pairing_3])

    f = h5open(filename_pairing, "w")
    write(f, "pairing", pairing)
    close(f)
end