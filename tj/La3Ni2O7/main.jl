# tJ bilayer Square lattice  with symmetry condition
# 控制粒子数 和 total spin
#
using ITensors
using ITensors.HDF5
#
function tune_dopping(N, dop_level)
    L_cell = 2 * dop_level[2]    
    L_emp = 2 * dop_level[1]
    N_cell = div(N, L_cell)
    @show N_cell
    if N % L_cell > 0
        println("N / L_cell must be an integer")
        sqrt(-1)
    else
        state = fill("a", N)
        L_spin = L_cell - L_emp
        @show L_spin
        if L_spin % 2 > 0
            print("L_spin must be an even to keep total Sz=0")
            sqrt(-1)
        else
            for it = 1:N_cell 
                for j = 1:L_spin
                    if j % 2 ==1
                        state[(it-1) * L_cell + j] = "Up"
                        #@show (it-1) * L_cell + j
                    else
                        state[(it-1) * L_cell + j] = "Dn"
                        #@show (it-1) * L_cell + j
                    end
                end
                for jj = L_spin+1 : L_cell
                    state[(it-1) * L_cell + jj] = "Emp"
                    #@show (it-1) * L_cell + jj
                end
            end
        end
        return state
    end
end
#
#  分情况写晶格的格点编号
function lattice_z(Ly, Lx) # along z direaction   
    b = fill(0, (Ly * Lx, 2)) # for Jz direaction
    for i =1:Ly*Lx
        b[i,1] = (2 * i - 1)
        b[i,2] = 2 * i
    end
    println("bonds along z direaction")
    @show b
    return b
end
#
function lattice_y_uplayer(Lz, Ly, Lx) # along y direaction up layer
    #
    b = fill(0, ((Ly-1) * Lx,2))
    for i = 1:Lx
        for j =1:Ly-1
            b[(i-1) * (Ly-1) + j, 1] = (i-1) * Lz * Ly + (j-1) * Lz + 1
            b[(i-1) * (Ly-1) + j, 2] = (i-1) * Lz * Ly + (j-1) * Lz + 1 + Lz
        end
    end
    println("bonds along y direaction, up layer")
    @show b
    return b 
end
#
function lattice_y_uplayer_boundary(Lz, Ly, Lx) # Y boundary on up layer
    #
    b = fill(0, (Lx,2))
    for i = 1:Lx
        b[i,1] = (i-1) * Lz * Ly + 1
        b[i,2] = (i-1) * Lz * Ly + 1 + (Ly-1) * Lz 
    end
    println("bonds Y pbc, up layer")
    @show b
    return b
end
#
function lattice_y_downlayer(Lz, Ly, Lx) # along y direaction down layer
    #
    b = fill(0, ((Ly-1)*Lx,2))
    for i = 1:Lx
        for j = 1:Ly-1
            b[(i-1) * (Ly-1) + j, 1] = (i-1) * Lz * Ly + (j-1) * Lz + 2
            b[(i-1) * (Ly-1) + j, 2] = (i-1) * Lz * Ly + (j-1) * Lz + 2 + Lz 
        end
    end
    println("bonds along y direaction, down layer")
    @show b
    return b 
end
#
function lattice_y_downlayer_boundary(Lz, Ly, Lx) # Y boundary on down layer
    #
    b = fill(0, (Lx,2))
    for i = 1:Lx
        b[i,1] = (i-1) * Lz * Ly + 1 + 1
        b[i,2] = (i-1) * Lz * Ly + 1 + 1 + (Ly-1) * Lz
    end
    println("bonds Y pbc, down layer")
    @show b
    return b
end
#
function lattice_x_uplayer(Lz, Ly, Lx) # along x direaction up layer
    #
    b = fill(0, (Ly * (Lx-1),2))
    for i = 1:Lx-1 
        for j = 1:Ly
            b[(i-1) * (Ly) + j, 1] = (i-1) * Lz * Ly + (j-1) * Lz + 1 
            b[(i-1) * (Ly) + j, 2] = (i-1) * Lz * Ly + (j-1) * Lz + 1 + Lz * Ly
        end
    end
    println("bonds along x direaction, up layer")
    @show b
    return b
end
#
function lattice_x_downlayer(Lz, Ly, Lx) # along x direaction down layer
    #
    b = fill(0, (Ly * (Lx-1),2))
    for i = 1:Lx-1
        for j = 1:Ly
            b[(i-1) * (Ly) + j, 1] = (i-1) * Lz * Ly + (j-1) * Lz + 2 
            b[(i-1) * (Ly) + j, 2] = (i-1) * Lz * Ly + (j-1) * Lz + 2 + Lz * Ly
        end
    end
    println("bonds along x direaction, down layer")
    @show b
    return b
end
#
function build_sites(N,Is_conserve_qns)   
    if Is_conserve_qns == true
        sites = siteinds("tJ",N; conserve_qns = true, conserve_nf=true, conserve_sz=true, conserve_nfparity=true) # 加上 conserve_qns 条件
    elseif Is_conserve_qns == false
        sites = siteinds("tJ", N)
    else
        println("wrong input of Is_conserve_qns")
        sqrt(-1)
    end
    return sites
end
#
function get_psi0(N, dop_level, sites, Is_conserve_qns,way_of_psi0,loadfile)
    if way_of_psi0 == "default"
        println("get psi0 by the way of default")
        # get the initial state psi0
        if Is_conserve_qns == true
            state = tune_dopping(N, dop_level)
            @show state
            println("initial state psi0 is random but with conserved sites")
            psi0 = randomMPS(sites, state, N)
            println("setting initial psi0")
            @show flux(psi0)
        elseif Is_conserve_qns == false
            state = tune_dopping(N, dop_level)
            @show state
            println("initial state psi0 is random !!!<without>!!! conserved sites")
            psi0 = randomMPS(sites, state, N)
            println("setting initial psi0")
            @show flux(psi0)
        else
            println("wrong input of Is_conserve_qns")
            sqrt(-1)
        end
        return psi0
    elseif way_of_psi0 == "readfile"
        println("get psi0 by the way of readfile")
        println("filename is:")
        println(loadfile)
        f = h5open(loadfile, "r")
        psi0 = read(f, "psi", MPS)
        close(f)
        readmaxdim = maxlinkdim(psi0)
        println("**********************************")
        println("read psi from disk, the maxdim is:")
        @show readmaxdim
        println("**********************************")
        return psi0
    else
        println("wrong input of way_of_psi0")
        sqrt(-1)
    end
end
#
function get_opsum(Yboundary, Lz, Ly, Lx,tz,Jz,Nz, txy,Jxy, Nxy)
    # get lattice bonds
    b_z = lattice_z(Ly, Lx)
    b_y1 = lattice_y_uplayer(Lz, Ly, Lx)
    b_y1_boundary = lattice_y_uplayer_boundary(Lz, Ly, Lx)
    b_y2 = lattice_y_downlayer(Lz, Ly, Lx)
    b_y2_boundary = lattice_y_downlayer_boundary(Lz, Ly, Lx)
    b_x1 = lattice_x_uplayer(Lz, Ly, Lx)
    b_x2 = lattice_x_downlayer(Lz, Ly, Lx)
    #
    os = OpSum()
    for j = 1:size(b_z)[1]
        os += -tz, "Cdagup", b_z[j,1], "Cup", b_z[j,2]
        os += -tz, "Cdagup", b_z[j,2], "Cup", b_z[j,1]

        os += -tz, "Cdagdn", b_z[j,1], "Cdn", b_z[j,2]
        os += -tz, "Cdagdn", b_z[j,2], "Cdn", b_z[j,1]

        os += Jz, "Sz", b_z[j,1], "Sz", b_z[j,2]
        os += 0.5 * Jz, "S+", b_z[j,1], "S-", b_z[j,2]
        os += 0.5 * Jz, "S-", b_z[j,1], "S+", b_z[j,2] 
        os += Nz, "Ntot", b_z[j,1], "Ntot", b_z[j,2]
    end
    for j = 1:size(b_y1)[1]
        os += -txy, "Cdagup", b_y1[j,1], "Cup", b_y1[j,2]
        os += -txy, "Cdagup", b_y1[j,2], "Cup", b_y1[j,1]

        os += -txy, "Cdagdn", b_y1[j,1], "Cdn", b_y1[j,2]
        os += -txy, "Cdagdn", b_y1[j,2], "Cdn", b_y1[j,1]

        os += Jxy, "Sz", b_y1[j,1], "Sz", b_y1[j,2]
        os += 0.5 * Jxy, "S+", b_y1[j,1], "S-", b_y1[j,2]
        os += 0.5 * Jxy, "S-", b_y1[j,1], "S+", b_y1[j,2] 
        os += Nxy, "Ntot", b_y1[j,1], "Ntot", b_y1[j,2]
    end

    for j = 1:size(b_y2)[1]
        os += -txy, "Cdagup", b_y2[j,1], "Cup", b_y2[j,2]
        os += -txy, "Cdagup", b_y2[j,2], "Cup", b_y2[j,1]

        os += -txy, "Cdagdn", b_y2[j,1], "Cdn", b_y2[j,2]
        os += -txy, "Cdagdn", b_y2[j,2], "Cdn", b_y2[j,1]

        os += Jxy, "Sz", b_y2[j,1], "Sz", b_y2[j,2]
        os += 0.5 * Jxy, "S+", b_y2[j,1], "S-", b_y2[j,2]
        os += 0.5 * Jxy, "S-", b_y2[j,1], "S+", b_y2[j,2] 
        os += Nxy, "Ntot", b_y2[j,1], "Ntot", b_y2[j,2]
    end

    for j = 1:size(b_x1)[1]
        os += -txy, "Cdagup", b_x1[j,1], "Cup", b_x1[j,2]
        os += -txy, "Cdagup", b_x1[j,2], "Cup", b_x1[j,1]

        os += -txy, "Cdagdn", b_x1[j,1], "Cdn", b_x1[j,2]
        os += -txy, "Cdagdn", b_x1[j,2], "Cdn", b_x1[j,1]

        os += Jxy, "Sz", b_x1[j,1], "Sz", b_x1[j,2]
        os += 0.5 * Jxy, "S+", b_x1[j,1], "S-", b_x1[j,2]
        os += 0.5 * Jxy, "S-", b_x1[j,1], "S+", b_x1[j,2] 
        os += Nxy, "Ntot", b_x1[j,1], "Ntot", b_x1[j,2]
    end
    for j = 1:size(b_x2)[1]
        os += -txy, "Cdagup", b_x2[j,1], "Cup", b_x2[j,2]
        os += -txy, "Cdagup", b_x2[j,2], "Cup", b_x2[j,1]

        os += -txy, "Cdagdn", b_x2[j,1], "Cdn", b_x2[j,2]
        os += -txy, "Cdagdn", b_x2[j,2], "Cdn", b_x2[j,1]

        os += Jxy, "Sz", b_x2[j,1], "Sz", b_x2[j,2]
        os += 0.5 * Jxy, "S+", b_x2[j,1], "S-", b_x2[j,2]
        os += 0.5 * Jxy, "S-", b_x2[j,1], "S+", b_x2[j,2] 
        os += Nxy, "Ntot", b_x2[j,1], "Ntot", b_x2[j,2]
    end
    if Yboundary == "pbc"
        println("Yboundary: Periodic boundary condition")
        for j = 1:size(b_y1_boundary)[1]
            os += -txy, "Cdagup", b_y1_boundary[j,1], "Cup", b_y1_boundary[j,2]
            os += -txy, "Cdagup", b_y1_boundary[j,2], "Cup", b_y1_boundary[j,1]

            os += -txy, "Cdagdn", b_y1_boundary[j,1], "Cdn", b_y1_boundary[j,2]
            os += -txy, "Cdagdn", b_y1_boundary[j,2], "Cdn", b_y1_boundary[j,1]

            os += Jxy, "Sz", b_y1_boundary[j,1], "Sz", b_y1_boundary[j,2]
            os += 0.5 * Jxy, "S+", b_y1_boundary[j,1], "S-", b_y1_boundary[j,2]
            os += 0.5 * Jxy, "S-", b_y1_boundary[j,1], "S+", b_y1_boundary[j,2] 
            os += Nxy, "Ntot", b_y1_boundary[j,1], "Ntot", b_y1_boundary[j,2]
        end
        for j = 1:size(b_y2_boundary)[1]
            os += -txy, "Cdagup", b_y2_boundary[j,1], "Cup", b_y2_boundary[j,2]
            os += -txy, "Cdagup", b_y2_boundary[j,2], "Cup", b_y2_boundary[j,1]

            os += -txy, "Cdagdn", b_y2_boundary[j,1], "Cdn", b_y2_boundary[j,2]
            os += -txy, "Cdagdn", b_y2_boundary[j,2], "Cdn", b_y2_boundary[j,1]

            os += Jxy, "Sz", b_y2_boundary[j,1], "Sz", b_y2_boundary[j,2]
            os += 0.5 * Jxy, "S+", b_y2_boundary[j,1], "S-", b_y2_boundary[j,2]
            os += 0.5 * Jxy, "S-", b_y2_boundary[j,1], "S+", b_y2_boundary[j,2] 
            os += Nxy, "Ntot", b_y2_boundary[j,1], "Ntot", b_y2_boundary[j,2]
        end
    elseif Yboundary == "obc"
        println("Yboundary: Open boundary condition")
    else
        println("wrong input of Yboundary")
        sqrt(-1)
    end
    return os
end
#
function rundmrg()
    #-----------------------------------------------------------------------
    # control parameter
    bc = "pbc" # {pbc, obc} along Y direaction
    Lz = 2
    Ly = 3
    Lx = 8
    N = Lz * Lx * Ly
    #
    txy = 3.0  # txy: t on xy plane
    Jxy = 1.0  # Jxy: J on xy plane
    Nxy = -0.25 * Jxy   # Nxy: for the n_in_j interaction

    tz = 0.0  # tz: t between two layers
    Jz = 0.5  # Jz: J between two layers
    Nz = 0.0 * Jz  # Nz: for the n_in_j interaction 
    dop_level = (1,2) # delta    
    #
    Is_conserve_qns = true
    way_of_psi0_list = ["default", "readfile"] # {readfile or default}
    nepoch = 3
    #
    nsweeps = 4
    maxdim_list = [[100, 150, 200,200],[200,300,300,400],[400,500,600,600]] # pre setting cutoff dimemsion
    write_disk_dim = 3000 # 
    cutoff = [1E-5]
    noise = [1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-11, 0.0]
    #
    workpath = "./output/psi/"
    #
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_txy$(txy)_tz$(tz)_Jxy$(Jxy)_Jz$(Jz)_dop$(dop_level)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    #------------------------------------------------------------------------------
    os = get_opsum(bc, Lz, Ly, Lx,tz,Jz,Nz, txy,Jxy, Nxy)
    #-------------------------------------------------------------------------------
    # begin dmrg
    for in = 1:nepoch
        println("**************begin $(in) sweeps epoch************************")
        if in == 1
            way_of_psi0 = way_of_psi0_list[1]
        else
            way_of_psi0 = way_of_psi0_list[2]
        end
        # get the hamiltonian MPO
        sites = build_sites(N, Is_conserve_qns)
        psi0 = get_psi0(N, dop_level,sites,Is_conserve_qns,way_of_psi0,filename_psi) # 
        sites = siteinds(psi0)    
        H = MPO(os, sites)
        #
        maxdim = maxdim_list[in]
        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise,write_when_maxdim_exceeds=write_disk_dim)
        #
        println("------------------show energy-----------------------------")
        @show energy
        @show flux(psi)
        @show maxlinkdim(psi)
        #
        @show (Lz, Ly, Lx, txy, Jxy, tz, Jz, dop_level)
        per_energy = energy / N
        @show per_energy
        #Compute the energy variance of an MPS to check whether it is an eigenstate.
        #H2 = inner(H,psi,H,psi)
        #E = inner(psi',H,psi)
        #@time var = H2-E^2
        #@show var
        println("------------------save psi-----------------------------")
        #psi0 = psi            
        #-----------------------------------------------------
        # save psi
        f = h5open(filename_psi, "w")
        write(f, "psi", psi)
        close(f)
        #end
    end
    return
end
#
let 
    rundmrg()
end
#
