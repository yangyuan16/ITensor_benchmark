#
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#
if __name__ == "__main__":
    print()
    bc = "pbc"
    Lz = 2
    Ly = 3
    Lx = 48
    N = Lz * Ly * Lx
    S = 0.5
    dop = 0.5
    t = 3
    J = 1
    Jz = 0.5
    dim = 4000 # dim cutoff
    workpath = "E:/WORK/Work/Project/La3Ni2O7"
    filepath1 = "/data_itensor/Lz%d_Ly%d_Lx%d/dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "/t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filename = workpath + filepath1 + filepath2 + "/corrz.h5"
    print("filename: ", filename)
    hdfFile = h5py.File(filename, 'r')
    data = hdfFile.get('corrz')
    print("HDF5 data: ", data)
    corrz = list(data)
    print(corrz)
    hdfFile.close()
    df = pd.DataFrame(columns=["site1", "site2", "corre"])
    df["site1"] = corrz[0]
    df["site2"] = corrz[1]
    df["corre"] = corrz[2]
    df = df[(df["site2"]-df["site1"]) % (Lz * Ly) == 0] 
    df.sort_values(["site2"],inplace=True)
    print(df.head())
    print(df.tail())
    print(len(df))
    #==============================================
    site1 = df["site1"].values
    site2 = df["site2"].values
    corre = np.abs(np.array(df["corre"].values))
    r = np.array((site2 - site1) / (Lz * Ly))
    print(r)
    #-----------------plot ax1-------------------------
    fig = plt.figure(figsize=(6,10)) 
    ax1 = plt.subplot(3,1,1)
    ax2 = plt.subplot(3,1,2)
    ax3 = plt.subplot(3,1,3)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    text = "Lz%d_Ly%d_Lx%d_dop%g_" % (Lz, Ly, Lx,dop) + "t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    ax1.text(0.3, 0.009, text)
    ax1.plot(r,corre,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"|i-j|"
    label_y = "Spin Correlation"
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    #-----------------plot ax2--------------------------
    plt.sca(ax2)  ##选择对ax2进行绘图
    ax2=plt.gca() #获得坐标轴的句柄
    #corre1 = np.log10(corre)
    ax2.plot(r,corre,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"|i-j|"
    label_y = "log(Spin Correlation)"
    plt.yscale("log")
    ax2.set_xlabel(label_x, size= 14)
    ax2.set_ylabel(label_y, size= 14)
    ax2.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax2.set_xlim([0,8])
    #ax2.set_ylim([-0.1,1])
    #ax2.set_xticks([0,2,4,6,8])
    #ax2.set_yticks([-0.1,0,0.5,1])
    #ax2.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    #----------------plot ax3---------------------------
    plt.sca(ax3)  ##选择对ax2进行绘图
    ax3=plt.gca() #获得坐标轴的句柄
    #corre1 = np.log10(corre)
    #r1 = np.log(r)
    ax3.plot(r,corre,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    plt.yscale("log")
    plt.xscale("log")
    label_x = r"log(|i-j|)"
    label_y = r"log(Spin Correlation)"
    ax3.set_xlabel(label_x, size= 14)
    ax3.set_ylabel(label_y, size= 14)
    ax3.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax3.set_xlim([0,8])
    #ax3.set_ylim([-0.1,1])
    #ax3.set_xticks([0,2,4,6,8])
    #ax3.set_yticks([-0.1,0,0.5,1])
    #ax3.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    #
    plt.show()
    







