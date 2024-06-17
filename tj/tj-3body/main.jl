#
using ITensors
using ITensors.HDF5
# 风情况写格点编号
function lattice_x(Ly,Lx) # along x direaction
    b = fill(0, ((Lx-1)*Ly,2)) # along x direaction
    for i = 1:Lx-1
          for j = 1:Ly
              b[(i-1) * (Ly) + j, 1] = (i-1) * Ly + j
              b[(i-1) * (Ly) + j, 2] = (i-1) * Ly + j + Ly
          end
    end
    println("bonds along x direaction")
    @show b
    return b
end
#
function lattice_y(Ly, Lx) # along y direaction
    b = fill(0, ((Ly-1)*Lx,2)) # along y direaction
    for i = 1:Lx
        for j = 1:Ly-1
          b[(i-1)*(Ly-1)+j, 1] = (i-1) * Ly + j
          b[(i-1)*(Ly-1)+j, 2] = (i-1) * Ly + j + 1           
        end
    end
    println("bonds along y direaction")
    @show b
    return b
end
#
function lattice_y_boundary(Ly, Lx) # Y boundary for two body interaction
    #
    b = fill(0, (Lx,2))
    for i = 1:Lx
        b[i,1] = (i-1) * Ly + 1
        b[i,2] = (i-1) * Ly + Ly
    end
    println("bonds Y pbc")
    @show b
    return b      
end
#
function lattice_ijk_1(Ly, Lx)
    #=
    i------ 
    |     |
    j-----k
    =#
    b = fill(0, ((Ly-1)*(Lx-1),3))
    for i = 1:(Lx-1)
          for j = 1:(Ly-1)
              #println((i-1) * (Ly) + j)
              b[(i-1) * (Ly-1) + j, 1] = (i-1) * Ly + j + 1    
              b[(i-1) * (Ly-1) + j, 2] = (i-1) * Ly + j
              b[(i-1) * (Ly-1) + j, 3] = (i-1) * Ly + j + Ly
          end
    end
    println("3body ijk bonds of type 1")
    @show b
    return b
end
#
function lattice_ijk_1_boundary(Ly, Lx)
          #=
          i------ 
          |     |
          j-----k
          =#
          b = fill(0, ((Lx-1),3))
          for i = 1:(Lx-1)
              b[i, 1] = (i-1) * Ly + Ly    
              b[i, 2] = (i-1) * Ly + 1
              b[i, 3] = (i-1) * Ly + 1 + Ly
          end
          println("3body ijk bonds PBC Boundary of type 1")
          @show b
          return b
end
#
function lattice_ijk_2(Ly, Lx)
          #=
          ------k 
          |     |
          i-----j
          =#
          b = fill(0, ((Ly-1)*(Lx-1),3))
          for i = 1:(Lx-1)
                for j = 1:(Ly-1)
                    #println((i-1) * (Ly) + j)
                    b[(i-1) * (Ly-1) + j, 1] = (i-1) * Ly + j    
                    b[(i-1) * (Ly-1) + j, 2] = (i-1) * Ly + j + Ly
                    b[(i-1) * (Ly-1) + j, 3] = (i-1) * Ly + j + Ly + 1
                end
          end
          println("3body ijk bonds of type 2")
          @show b
          return b
end
#
function lattice_ijk_2_boundary(Ly, Lx)
          #=
          ------k 
          |     |
          i-----j
          =#
          b = fill(0, ((Lx-1),3))
          for i = 1:(Lx-1)
              b[i, 1] = (i-1) * Ly + 1   
              b[i, 2] = (i-1) * Ly + 1 + Ly
              b[i, 3] = (i-1) * Ly + Ly + Ly
          end
          println("3body ijk bonds PBC Boundary of type 2")
          @show b
          return b
end
#
function lattice_ijk_3(Ly, Lx)
          #=
          k-----j 
          |     |
          ------i
          =#
          b = fill(0, ((Ly-1)*(Lx-1),3))
          for i = 1:(Lx-1)
                for j = 1:(Ly-1)
                    #println((i-1) * (Ly) + j)
                    b[(i-1) * (Ly-1) + j, 1] = (i-1) * Ly + j + Ly   
                    b[(i-1) * (Ly-1) + j, 2] = (i-1) * Ly + j + Ly + 1
                    b[(i-1) * (Ly-1) + j, 3] = (i-1) * Ly + j + 1
                end
          end
          println("3body ijk bonds of type 3")
          @show b
          return b
end
#
function lattice_ijk_3_bounday(Ly, Lx)
          #=
          k-----j 
          |     |
          ------i
          =#
          b = fill(0, ((Lx-1),3))
          for i = 1:(Lx-1)
              #println((i-1) * (Ly) + j)
              b[i, 1] = (i-1) * Ly + Ly + 1  
              b[i, 2] = (i-1) * Ly + Ly + Ly
              b[i, 3] = (i-1) * Ly + Ly
          end
          println("3body ijk bonds PBC Boundary of type 3")
          @show b
          return b
end
#
function lattice_ijk_4(Ly, Lx)
          #=
          j-----i 
          |     |
          k------
          =#
          b = fill(0, ((Ly-1)*(Lx-1),3))
          for i = 1:(Lx-1)
                for j = 1:(Ly-1)
                    #println((i-1) * (Ly) + j)
                    b[(i-1) * (Ly-1) + j, 1] = (i-1) * Ly + j + Ly + 1 
                    b[(i-1) * (Ly-1) + j, 2] = (i-1) * Ly + j + 1
                    b[(i-1) * (Ly-1) + j, 3] = (i-1) * Ly + j 
                end
          end
          println("3body ijk bonds of type 4")
          @show b
          return b
end
#
function lattice_ijk_4_boundary(Ly, Lx)
          #=
          j-----i 
          |     |
          k------
          =#
          b = fill(0, (Lx-1,3))
          for i = 1:(Lx-1)
              #println((i-1) * (Ly) + j)
              b[i, 1] = (i-1) * Ly + Ly + Ly 
              b[i, 2] = (i-1) * Ly + Ly
              b[i, 3] = (i-1) * Ly + 1 
          end
          println("3body ijk bonds PBC Boundary of type 4")
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
function get_state(N,doph)
    N_group = Int(floor(N / doph))
    N_res = N % doph
    state = fill("a",N)
    if (N-doph) % 2 ==1
       println("N-doph is a odd number, can not statisfy Total Sz = 0 condition")
       sqrt(-1)
    elseif N_group > 1   
        for i =1:doph
            for j = 1:N_group
                    k = (i-1) * N_group + j - 1
                if k % N_group == 0
                    state[k+1] = "Emp"
                elseif (i % 2 == 1) && ((k % N_group)%2 ==1)
                    state[k+1] = "Up"
                elseif (i % 2 == 1) && ((k % N_group)%2 ==0)
                    state[k+1] = "Dn"
                elseif (i % 2 == 0) && ((k % N_group)%2 ==1)
                    state[k+1] = "Dn"
                elseif (i % 2 == 0) && ((k % N_group)%2 ==0)
                    state[k+1] = "Up"
                end
            end
        end
        for i = 0:N_res-1
            if ((doph) % 2 == 0) && ((i % N_group)%2==1)
                state[doph*N_group+i+1] = "Dn"
            elseif ((doph) % 2 == 0) && ((i % N_group)%2==0)
                state[doph*N_group+i+1] = "Up"
            elseif ((doph) % 2 == 1) && ((i % N_group)%2==1)
                state[doph*N_group+i+1] = "Up"
            elseif ((doph) % 2 == 1) && ((i % N_group)%2==0)
                state[doph*N_group+i+1] = "Dn"            
            end
        end
    elseif N_group == 1
        N_L = Int((N - doph) / 2)  
        for i = 1:N_L
            state[i] = "Up"
        end
        for i = (N_L+1):(N_L+doph)
            state[i] = "Emp"
        end
        for i = (N_L+doph+1):N
            state[i] = "Dn"
        end          
    end
    return state
end
#


#


#


#


#





#
#
function get_opsum(Yboundary, Ly, Lx, t, J, Nij, Jijk)
    # get lattice bonds
    b_x = lattice_x(Ly,Lx)
    b_y = lattice_y(Ly,Lx)
    b_y_pbc = lattice_y_boundary(Ly,Lx)
    b_ijk_1 = lattice_ijk_1(Ly,Lx)
    b_ijk_1_pbc = lattice_ijk_1_boundary(Ly,Lx)
    b_ijk_2 = lattice_ijk_2(Ly,Lx)
    b_ijk_2_pbc = lattice_ijk_2_boundary(Ly,Lx)
    b_ijk_3 = lattice_ijk_3(Ly,Lx)
    b_ijk_3_pbc = lattice_ijk_3_bounday(Ly,Lx)
    b_ijk_4 = lattice_ijk_4(Ly,Lx)
    b_ijk_4_pbc = lattice_ijk_4_boundary(Ly,Lx)

    os = OpSum()
    for j = 1:size(b_x)[1]
        os += -t, "Cdagup", b_x[j,1], "Cup", b_x[j,2]
        os += -t, "Cdagup", b_x[j,2], "Cup", b_x[j,1]

        os += -t, "Cdagdn", b_x[j,1], "Cdn", b_x[j,2]
        os += -t, "Cdagdn", b_x[j,2], "Cdn", b_x[j,1]

        os += J, "Sz", b_x[j,1], "Sz", b_x[j,2]
        os += 0.5 * J, "S+", b_x[j,1], "S-", b_x[j,2]
        os += 0.5 * J, "S-", b_x[j,1], "S+", b_x[j,2]
        os += Nij, "Ntot", b_x[j,1], "Ntot", b_x[j,2]
    end
#
    for j = 1:size(b_y)[1]
          os += -t, "Cdagup", b_y[j,1], "Cup", b_y[j,2]
          os += -t, "Cdagup", b_y[j,2], "Cup", b_y[j,1]
  
          os += -t, "Cdagdn", b_y[j,1], "Cdn", b_y[j,2]
          os += -t, "Cdagdn", b_y[j,2], "Cdn", b_y[j,1]
  
          os += J, "Sz", b_y[j,1], "Sz", b_y[j,2]
          os += 0.5 * J, "S+", b_y[j,1], "S-", b_y[j,2]
          os += 0.5 * J, "S-", b_y[j,1], "S+", b_y[j,2]
          os += Nij, "Ntot", b_y[j,1], "Ntot", b_y[j,2]
      end
#
    for j = 1:size(b_ijk_1)[1]
          os += Jijk, "Cdagup", b_ijk_1[j,1], "Cdagdn", b_ijk_1[j,2], "Cdn", b_ijk_1[j,2], "Cup", b_ijk_1[j,3]
          os += Jijk, "Cdagup", b_ijk_1[j,3], "Cdagdn", b_ijk_1[j,2], "Cdn", b_ijk_1[j,2], "Cup", b_ijk_1[j,1]

          os += Jijk, "Cdagdn", b_ijk_1[j,1], "Cdagup", b_ijk_1[j,2], "Cup", b_ijk_1[j,2], "Cdn", b_ijk_1[j,3]
          os += Jijk, "Cdagdn", b_ijk_1[j,3], "Cdagup", b_ijk_1[j,2], "Cup", b_ijk_1[j,2], "Cdn", b_ijk_1[j,1]

          os += -Jijk, "Cdagup", b_ijk_1[j,1], "Cdagdn", b_ijk_1[j,2], "Cup", b_ijk_1[j,2], "Cdn", b_ijk_1[j,3]
          os += -Jijk, "Cdagup", b_ijk_1[j,3], "Cdagdn", b_ijk_1[j,2], "Cup", b_ijk_1[j,2], "Cdn", b_ijk_1[j,1]

          os += -Jijk, "Cdagdn", b_ijk_1[j,1], "Cdagup", b_ijk_1[j,2], "Cdn", b_ijk_1[j,2], "Cup", b_ijk_1[j,3]
          os += -Jijk, "Cdagdn", b_ijk_1[j,3], "Cdagup", b_ijk_1[j,2], "Cdn", b_ijk_1[j,2], "Cup", b_ijk_1[j,1]
    end 
#
    for j = 1:size(b_ijk_2)[1]
          os += Jijk, "Cdagup", b_ijk_2[j,1], "Cdagdn", b_ijk_2[j,2], "Cdn", b_ijk_2[j,2], "Cup", b_ijk_2[j,3]
          os += Jijk, "Cdagup", b_ijk_2[j,3], "Cdagdn", b_ijk_2[j,2], "Cdn", b_ijk_2[j,2], "Cup", b_ijk_2[j,1]

          os += Jijk, "Cdagdn", b_ijk_2[j,1], "Cdagup", b_ijk_2[j,2], "Cup", b_ijk_2[j,2], "Cdn", b_ijk_2[j,3]
          os += Jijk, "Cdagdn", b_ijk_2[j,3], "Cdagup", b_ijk_2[j,2], "Cup", b_ijk_2[j,2], "Cdn", b_ijk_2[j,1]

          os += -Jijk, "Cdagup", b_ijk_2[j,1], "Cdagdn", b_ijk_2[j,2], "Cup", b_ijk_2[j,2], "Cdn", b_ijk_2[j,3]
          os += -Jijk, "Cdagup", b_ijk_2[j,3], "Cdagdn", b_ijk_2[j,2], "Cup", b_ijk_2[j,2], "Cdn", b_ijk_2[j,1]

          os += -Jijk, "Cdagdn", b_ijk_2[j,1], "Cdagup", b_ijk_2[j,2], "Cdn", b_ijk_2[j,2], "Cup", b_ijk_2[j,3]
          os += -Jijk, "Cdagdn", b_ijk_2[j,3], "Cdagup", b_ijk_2[j,2], "Cdn", b_ijk_2[j,2], "Cup", b_ijk_2[j,1]
    end 
#
    for j = 1:size(b_ijk_3)[1]
          os += Jijk, "Cdagup", b_ijk_3[j,1], "Cdagdn", b_ijk_3[j,2], "Cdn", b_ijk_3[j,2], "Cup", b_ijk_3[j,3]
          os += Jijk, "Cdagup", b_ijk_3[j,3], "Cdagdn", b_ijk_3[j,2], "Cdn", b_ijk_3[j,2], "Cup", b_ijk_3[j,1]

          os += Jijk, "Cdagdn", b_ijk_3[j,1], "Cdagup", b_ijk_3[j,2], "Cup", b_ijk_3[j,2], "Cdn", b_ijk_3[j,3]
          os += Jijk, "Cdagdn", b_ijk_3[j,3], "Cdagup", b_ijk_3[j,2], "Cup", b_ijk_3[j,2], "Cdn", b_ijk_3[j,1]

          os += -Jijk, "Cdagup", b_ijk_3[j,1], "Cdagdn", b_ijk_3[j,2], "Cup", b_ijk_3[j,2], "Cdn", b_ijk_3[j,3]
          os += -Jijk, "Cdagup", b_ijk_3[j,3], "Cdagdn", b_ijk_3[j,2], "Cup", b_ijk_3[j,2], "Cdn", b_ijk_3[j,1]

          os += -Jijk, "Cdagdn", b_ijk_3[j,1], "Cdagup", b_ijk_3[j,2], "Cdn", b_ijk_3[j,2], "Cup", b_ijk_3[j,3]
          os += -Jijk, "Cdagdn", b_ijk_3[j,3], "Cdagup", b_ijk_3[j,2], "Cdn", b_ijk_3[j,2], "Cup", b_ijk_3[j,1]
    end
#
    for j = 1:size(b_ijk_4)[1]
          os += Jijk, "Cdagup", b_ijk_4[j,1], "Cdagdn", b_ijk_4[j,2], "Cdn", b_ijk_4[j,2], "Cup", b_ijk_4[j,3]
          os += Jijk, "Cdagup", b_ijk_4[j,3], "Cdagdn", b_ijk_4[j,2], "Cdn", b_ijk_4[j,2], "Cup", b_ijk_4[j,1]

          os += Jijk, "Cdagdn", b_ijk_4[j,1], "Cdagup", b_ijk_4[j,2], "Cup", b_ijk_4[j,2], "Cdn", b_ijk_4[j,3]
          os += Jijk, "Cdagdn", b_ijk_4[j,3], "Cdagup", b_ijk_4[j,2], "Cup", b_ijk_4[j,2], "Cdn", b_ijk_4[j,1]

          os += -Jijk, "Cdagup", b_ijk_4[j,1], "Cdagdn", b_ijk_4[j,2], "Cup", b_ijk_4[j,2], "Cdn", b_ijk_4[j,3]
          os += -Jijk, "Cdagup", b_ijk_4[j,3], "Cdagdn", b_ijk_4[j,2], "Cup", b_ijk_4[j,2], "Cdn", b_ijk_4[j,1]

          os += -Jijk, "Cdagdn", b_ijk_4[j,1], "Cdagup", b_ijk_4[j,2], "Cdn", b_ijk_4[j,2], "Cup", b_ijk_4[j,3]
          os += -Jijk, "Cdagdn", b_ijk_4[j,3], "Cdagup", b_ijk_4[j,2], "Cdn", b_ijk_4[j,2], "Cup", b_ijk_4[j,1]
    end
#
    if Yboundary == "pbc"
          println("Yboundary: Periodic boundary condition")
          for j = 1:size(b_y_pbc)[1]
              os += -t, "Cdagup", b_y_pbc[j,1], "Cup", b_y_pbc[j,2]
              os += -t, "Cdagup", b_y_pbc[j,2], "Cup", b_y_pbc[j,1]
            
              os += -t, "Cdagdn", b_y_pbc[j,1], "Cdn", b_y_pbc[j,2]
              os += -t, "Cdagdn", b_y_pbc[j,2], "Cdn", b_y_pbc[j,1]

              os += J, "Sz", b_y_pbc[j,1], "Sz", b_y_pbc[j,2]
              os += 0.5 * J, "S+", b_y_pbc[j,1], "S-", b_y_pbc[j,2]
              os += 0.5 * J, "S-", b_y_pbc[j,1], "S+", b_y_pbc[j,2]
              os += Nij, "Ntot", b_y_pbc[j,1], "Ntot", b_y_pbc[j,2]
          end
          for j = 1:size(b_ijk_1_pbc)[1]
                    os += Jijk, "Cdagup", b_ijk_1_pbc[j,1], "Cdagdn", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,3]
                    os += Jijk, "Cdagup", b_ijk_1_pbc[j,3], "Cdagdn", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,1]
          
                    os += Jijk, "Cdagdn", b_ijk_1_pbc[j,1], "Cdagup", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,3]
                    os += Jijk, "Cdagdn", b_ijk_1_pbc[j,3], "Cdagup", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,1]
          
                    os += -Jijk, "Cdagup", b_ijk_1_pbc[j,1], "Cdagdn", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,3]
                    os += -Jijk, "Cdagup", b_ijk_1_pbc[j,3], "Cdagdn", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,1]
          
                    os += -Jijk, "Cdagdn", b_ijk_1_pbc[j,1], "Cdagup", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,3]
                    os += -Jijk, "Cdagdn", b_ijk_1_pbc[j,3], "Cdagup", b_ijk_1_pbc[j,2], "Cdn", b_ijk_1_pbc[j,2], "Cup", b_ijk_1_pbc[j,1]
              end 
          #
          for j = 1:size(b_ijk_2_pbc)[1]
                    os += Jijk, "Cdagup", b_ijk_2_pbc[j,1], "Cdagdn", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,3]
                    os += Jijk, "Cdagup", b_ijk_2_pbc[j,3], "Cdagdn", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,1]
          
                    os += Jijk, "Cdagdn", b_ijk_2_pbc[j,1], "Cdagup", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,3]
                    os += Jijk, "Cdagdn", b_ijk_2_pbc[j,3], "Cdagup", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,1]
          
                    os += -Jijk, "Cdagup", b_ijk_2_pbc[j,1], "Cdagdn", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,3]
                    os += -Jijk, "Cdagup", b_ijk_2_pbc[j,3], "Cdagdn", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,1]
          
                    os += -Jijk, "Cdagdn", b_ijk_2_pbc[j,1], "Cdagup", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,3]
                    os += -Jijk, "Cdagdn", b_ijk_2_pbc[j,3], "Cdagup", b_ijk_2_pbc[j,2], "Cdn", b_ijk_2_pbc[j,2], "Cup", b_ijk_2_pbc[j,1]
          end 
          #
          for j = 1:size(b_ijk_3_pbc)[1]
                    os += Jijk, "Cdagup", b_ijk_3_pbc[j,1], "Cdagdn", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,3]
                    os += Jijk, "Cdagup", b_ijk_3_pbc[j,3], "Cdagdn", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,1]
          
                    os += Jijk, "Cdagdn", b_ijk_3_pbc[j,1], "Cdagup", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,3]
                    os += Jijk, "Cdagdn", b_ijk_3_pbc[j,3], "Cdagup", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,1]
          
                    os += -Jijk, "Cdagup", b_ijk_3_pbc[j,1], "Cdagdn", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,3]
                    os += -Jijk, "Cdagup", b_ijk_3_pbc[j,3], "Cdagdn", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,1]
          
                    os += -Jijk, "Cdagdn", b_ijk_3_pbc[j,1], "Cdagup", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,3]
                    os += -Jijk, "Cdagdn", b_ijk_3_pbc[j,3], "Cdagup", b_ijk_3_pbc[j,2], "Cdn", b_ijk_3_pbc[j,2], "Cup", b_ijk_3_pbc[j,1]
          end
          #
          for j = 1:size(b_ijk_4_pbc)[1]
                    os += Jijk, "Cdagup", b_ijk_4_pbc[j,1], "Cdagdn", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,3]
                    os += Jijk, "Cdagup", b_ijk_4_pbc[j,3], "Cdagdn", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,1]
          
                    os += Jijk, "Cdagdn", b_ijk_4_pbc[j,1], "Cdagup", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,3]
                    os += Jijk, "Cdagdn", b_ijk_4_pbc[j,3], "Cdagup", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,1]
          
                    os += -Jijk, "Cdagup", b_ijk_4_pbc[j,1], "Cdagdn", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,3]
                    os += -Jijk, "Cdagup", b_ijk_4_pbc[j,3], "Cdagdn", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,1]
          
                    os += -Jijk, "Cdagdn", b_ijk_4_pbc[j,1], "Cdagup", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,3]
                    os += -Jijk, "Cdagdn", b_ijk_4_pbc[j,3], "Cdagup", b_ijk_4_pbc[j,2], "Cdn", b_ijk_4_pbc[j,2], "Cup", b_ijk_4_pbc[j,1]
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
#
function rundmrg()
    #-------------------------------------
    # control parameter
    bc = "pbc" # {pbc, obc} along Y direaction
    Ly = 4
    Lx = 8
    N = Lx * Ly
    
    t = 3.0 
    J = 1.0
    Nij = -0.25 * J # Nij: for the n_in_j interaction
    Jijk = -0.25 * J # Jijk: for three body interaction 
    #
    doph = 16 # number of holes
    #
    Is_conserve_qns = true
    nsweeps = 4
    maxdim = [100,200,300,400]
    cutoff = [1E-5]
    noise = [1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-11, 0.0]
    #
    workpath = "./"
    #
    filename_psi_1 = "Ly$(Ly)_Lx$(Lx)_N$(N)_t$(t)_J$(J)_dop$(doph)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    #------------------------------------------------------------------------------
    # begin dmrg
    # get the hamiltonian MPO
    sites = build_sites(N, Is_conserve_qns)
    os = get_opsum(bc, Ly, Lx, t, J, Nij, Jijk)
    state = get_state(N,doph)
    psi0 = randomMPS(sites, state, N)
    @show flux(psi0)
    #
    H = MPO(os, sites)
    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise)
    println("------------------show energy-----------------------------")
    @show energy
    @show flux(psi)
    @show maxlinkdim(psi)
    #
    @show (Ly, Lx, t, J, Jijk, doph)
    per_energy = energy / N
    @show per_energy
    #-----------------------------------------------------
    # save psi
    f = h5open(filename_psi, "w")
    write(f, "psi", psi)
    close(f)
end
#
let 
   rundmrg()
end