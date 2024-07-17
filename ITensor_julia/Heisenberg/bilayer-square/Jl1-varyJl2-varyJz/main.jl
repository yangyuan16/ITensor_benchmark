# Lz layerd squared lattice Heisenberg model
#
using ITensors
using ITensors.HDF5
# 
function lattice_obc(Lz, Ly, Lx)    
    b1 = fill(0, (Ly * Lx, 2)) # for Jz direaction
    for i =1:Ly*Lx
        b1[i,1] = (2 * i - 1)
        b1[i,2] = 2 * i
    end
    @show b1
    #
    b2 = fill(0, ((Lz*(Ly-1))*Lx, 2)) # Jy direaction for both up layer and down layer
    for i = 1:Lx
        s0 = (i-1) * Lz * Ly + 1
        for j = 1:(Lz * (Ly - 1))
            b2[(i-1) * (Lz * (Ly-1)) + j, 1] = s0 + (j-1)
            b2[(i-1) * (Lz * (Ly-1)) + j, 2] = s0 + (j+1)
        end
    end
    @show b2
    #
    b3 = fill(0, ((Lz * Ly) * (Lx-1),2)) # Jx direaction for both up layer and down layer
    for i = 1:(Lz*Ly)*(Lx-1)
        b3[i, 1] = i 
        b3[i, 2] = i + (Lz * Ly)
    end
    @show b3
    
    b12 = cat(b1,b2,dims=1)
    b123 = cat(b12,b3, dims=1)
    
    @show b123

    return b123
end 
#
function lattice_pbc(Lz, Ly, Lx)
    b1 = fill(0, (Ly * Lx, 2))
    for i =1:Ly*Lx
        b1[i,1] = (2 * i - 1)
        b1[i,2] = 2 * i
    end
    @show b1
    #
    b2 = fill(0, ((Lz*(Ly-1))*Lx, 2))
    for i = 1:Lx
        s0 = (i-1) * Lz * Ly + 1
        for j = 1:(Lz * (Ly - 1))
            b2[(i-1) * (Lz * (Ly-1)) + j, 1] = s0 + (j-1)
            b2[(i-1) * (Lz * (Ly-1)) + j, 2] = s0 + (j+1)
        end
    end
    @show b2
    #
    b3 = fill(0, ((Lz * Ly) * (Lx-1),2))
    for i = 1:(Lz*Ly)*(Lx-1)
        b3[i, 1] = i 
        b3[i, 2] = i + (Lz * Ly)
    end
    @show b3
    #
    b4 = fill(0, (Lx * Lz, 2))  # consider the periodic boundary
    for i = 1:Lx
        s0 = (i - 1)*(Lz * Ly) + 1
        for j = 1:Lz
            b4[(i-1) * Lz + j, 1] = s0 + j - 1
            b4[(i-1) * Lz + j, 2] = s0 + j - 1 + Lz * (Ly - 1) 
        end 
    end
    @show b4

    b12 = cat(b1,b2,dims=1)
    b123 = cat(b12,b3, dims=1)
    b1234 = cat(b123, b4, dims=1)

    return b1234

end

#  分情况写晶格的格点编号
function lattice_obc_z(Lz, Ly, Lx)    
      b1 = fill(0, (Ly * Lx, 2)) # for Jz direaction
      for i =1:Ly*Lx
          b1[i,1] = (2 * i - 1)
          b1[i,2] = 2 * i
      end
      @show b1
      return b1
end

function lattice_obc_y(Lz, Ly, Lx)
      #
      b2 = fill(0, ((Lz*(Ly-1))*Lx, 2)) # Jy direaction for both up layer and down layer
      for i = 1:Lx
          s0 = (i-1) * Lz * Ly + 1
          for j = 1:(Lz * (Ly - 1))
              b2[(i-1) * (Lz * (Ly-1)) + j, 1] = s0 + (j-1)
              b2[(i-1) * (Lz * (Ly-1)) + j, 2] = s0 + (j+1)
          end
      end
      @show b2
      return b2
end

function lattice_obc_x(Lz, Ly, Lx)
      #
      b3 = fill(0, ((Lz * Ly) * (Lx-1),2)) # Jx direaction for both up layer and down layer
      for i = 1:(Lz*Ly)*(Lx-1)
            b3[i, 1] = i 
            b3[i, 2] = i + (Lz * Ly)
      end
      @show b3
      return b3      
end

#
function rundmrg() 
    bc = "obc" # ob: open boundary condition; pb: periodic boundary condition
    Lz = 2
    Ly = 3
    Lx = 4
    N = Lz * Ly * Lx
    
    Jl1 = 1 #up layer
    Jl2 = 0.8 #down layer
    Jz = 0.6 # between two layers

    sites = siteinds("S=1/2", N; conserve_qns=true) # S=1/2 spins

    bz = lattice_obc_z(Lz, Ly, Lx)
    by = lattice_obc_y(Lz, Ly, Lx)
    bx = lattice_obc_x(Lz, Ly, Lx)

    
    os = OpSum()
    for j = 1:size(bz)[1]
        os += Jz, "Sz", bz[j,1], "Sz", bz[j,2]
        os += Jz * 0.5, "S+", bz[j,1], "S-", bz[j,2]
        os += Jz * 0.5, "S-", bz[j,1], "S+", bz[j,2]
    end

    for j = 1:size(by)[1]
      if by[j,1] % 2 == 1
            os += Jl1, "Sz", by[j,1], "Sz", by[j,2]
            os += Jl1 * 0.5, "S+", by[j,1], "S-", by[j,2]
            os += Jl1 * 0.5, "S-", by[j,1], "S+", by[j,2]
      else
            os += Jl2, "Sz", by[j,1], "Sz", by[j,2]
            os += Jl2 * 0.5, "S+", by[j,1], "S-", by[j,2]
            os += Jl2 * 0.5, "S-", by[j,1], "S+", by[j,2]
      end
    end
    
    for j = 1:size(bx)[1]
      if bx[j,1] % 2 == 1
            os += Jl1, "Sz", bx[j,1], "Sz", bx[j,2]
            os += Jl1 * 0.5, "S+", bx[j,1], "S-", bx[j,2]
            os += Jl1 * 0.5, "S-", bx[j,1], "S+", bx[j,2]
      else
            os += Jl2, "Sz", bx[j,1], "Sz", bx[j,2]
            os += Jl2 * 0.5, "S+", bx[j,1], "S-", bx[j,2]
            os += Jl2 * 0.5, "S-", bx[j,1], "S+", bx[j,2]
      end
    end


    
    H = MPO(os, sites)

    state = [isodd(n) ? "Up" : "Dn" for n = 1:N]
    psi0 = randomMPS(sites, state, 20)
    @show flux(psi0)

    nsweeps = 5
    maxdim = [10,20,100,100,200]
    cutoff = [1E-10]
    noise = [1E-6, 1E-7, 1E-8, 0.0]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)

    @show energy
    @show flux(psi)
    @show maxlinkdim(psi)
    per_energy = energy / N
    @show per_energy
    
    workpath = "./bilayer-square/Jl1-Jl2-varyJz/output/psi/"
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    f = h5open(filename_psi, "w")
    write(f, "psi", psi)
    close(f)
    
    return
end

let 
    @time rundmrg()
end