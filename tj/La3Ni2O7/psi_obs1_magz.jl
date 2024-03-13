using ITensors
using ITensors.HDF5
#
function magz(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
    workpath = "./output/psi/"
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_txy$(txy)_tz$(tz)_Jxy$(Jxy)_Jz$(Jz)_dop$(dop_level)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    #
    f = h5open(filename_psi, "r")
    psi = read(f, "psi", MPS)
    close(f)
    #--- measure magz ---
    magz = expect(psi, "Sz")
    @show magz
    savepath = "./"
    filename_magz_1 = filename_psi_1
    filename_magz_2 = filename_psi_2
    filename_magz_3 = "_magz.h5"
    filename_magz = join([savepath, filename_magz_1, filename_magz_2, filename_magz_3])

    f = h5open(filename_magz, "w")
    write(f, "magz", magz)
    close(f)
end
#
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
    println("---psi_obs1_magz---")
    magz(bc, Lz, Ly, Lx, N, txy, tz, Jxy, Jz, dop_level)
end
=#