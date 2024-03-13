import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
#
def plot_density_curve(r,density, label_y):
    fig = plt.figure(figsize=(10,3)) 
    ax1 = plt.subplot(1,1,1)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    text = "Lz%d_Ly%d_Lx%d_dop%g_" % (Lz, Ly, Lx,dop) + "t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    ax1.text(0.3, 0.009, text)
    ax1.plot(r,density,label=r"G",ls="-",lw=1.5,color="red",
             marker='o',alpha=1,markersize=6,markeredgewidth=1.5, markeredgecolor="k",
             markerfacecolor='w')
    label_x = r"i"
    label_y = label_y
    ax1.set_xlabel(label_x, size= 14)
    ax1.set_ylabel(label_y, size= 14)
    #
    ax1.tick_params(labelsize = 15) # 设置坐标刻度对应数字的大小
    #ax1.set_xlim([0,8])
    #ax1.set_ylim([-0.1,1])
    #ax1.set_xticks([0,2,4,6,8])
    #ax1.set_yticks([-0.1,0,0.5,1])
    #ax1.text(0.3,-0.05, r'$\mathrm{(a)}$', fontsize=18)
    plt.show()
    return
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
    filename = workpath + filepath1 + filepath2 + "/magz.h5"
    print("filename: ", filename)
    hdfFile = h5py.File(filename, 'r')
    data = hdfFile.get('magz')
    print("HDF5 data: ", data)
    magz = list(data)
    print(magz)
    hdfFile.close()
    #
    df = pd.DataFrame(columns=["site","magz"])
    df["site"] = range(1,N+1)
    df["magz"] = magz
    df.sort_values(['site'],inplace=True)
    df_layer1 = df[df["site"] % 2 == 0] 
    df_layer2 = df[df["site"] % 2 == 1]
    sites_layer1 = df_layer1["site"].values.reshape(-1,Ly).T
    magz_layer1 = df_layer1["magz"].values.reshape(-1, Ly).T
    sites_layer2 = df_layer2["site"].values.reshape(-1, Ly).T
    magz_layer2 = df_layer2["magz"].values.reshape(-1, Ly).T
    print(df.head())
    #
    #------------------------------------------------------------
    df_y0 = df[df['site'] % 6 ==0]
    r_y0 = df_y0["site"].values
    magz_y0 = df_y0["magz"].values
    
    df_y1 = df[df['site'] % 6 ==1]
    r_y1 = df_y1["site"].values
    magz_y1 = df_y1["magz"].values
    
    df_y2 = df[df['site'] % 6 ==2]
    r_y2 = df_y2["site"].values
    magz_y2 = df_y2["magz"].values

    df_y3 = df[df['site'] % 6 ==3]
    r_y3 = df_y3["site"].values
    magz_y3 = df_y3["magz"].values

    df_y4 = df[df['site'] % 6 ==4]
    r_y4 = df_y4["site"].values
    magz_y4 = df_y4["magz"].values

    df_y5 = df[df['site'] % 6 ==5]
    r_y5 = df_y5["site"].values
    magz_y5 = df_y5["magz"].values

    r_mean = range(len(r_y0))
    magz_mean = (np.array(magz_y0) + np.array(magz_y1) + np.array(magz_y2)  + 
                    np.array(magz_y3) + np.array(magz_y4) + np.array(magz_y5)) / 6
    #----------------------------------------------------------
    # plot the data
    fig = plt.figure(figsize=(12,2))
    ax1 = plt.subplot(2,1,1)
    ax2 = plt.subplot(2,1,2)
    #
    vmin1 = np.array(magz_layer1).min()
    vmax1 = np.array(magz_layer1).max()
    vmin2 = np.array(magz_layer2).min()
    vmax2 = np.array(magz_layer2).max()
    vmin = min(vmin1, vmin2)
    vmax = max(vmax1, vmax2)
    plt.sca(ax1)  ##选择对ax1进行绘图
    ax1=plt.gca() #获得坐标轴的句柄
    text = "Lz%d_Ly%d_Lx%d_dop%g_" % (Lz, Ly, Lx,dop) + "t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    plt.text(0.3, -0.25, text)
    sns.heatmap(magz_layer1, cmap="coolwarm", annot=True, linewidths=0.5, vmax=vmax, vmin=vmin,
                annot_kws={'size':9,'weight':'bold', 'color':'white', "rotation":90})
    ax1.set_title("magz <Si>, layer1")
    ax1.set_xlabel("")
    ax1.set_xticklabels([]) # 设置 x 轴图例为空置
    ax1.set_ylabel("")
    #---
    plt.sca(ax2)  ##选择对ax1进行绘图
    ax2=plt.gca() #获得坐标轴的句柄
    sns.heatmap(magz_layer2, cmap="coolwarm", annot=True, linewidths=0.5, vmax=vmax, vmin=vmin,
                annot_kws={'size':9,'weight':'bold', 'color':'white', "rotation":90})
    ax2.set_title("magz <Si>, layer2")
    ax2.set_xlabel("")
    ax2.set_xticklabels([]) # 设置 x 轴图例为空置
    ax2.set_ylabel("")
    #f.savefig("charge_density.jpg", bbox_inches="tight")
    plt.show()
    #----plot density curve----------------
    plot_density_curve(r=r_y0,density=magz_y0,label_y='<Si>(y0)')
    plot_density_curve(r=r_y1,density=magz_y1,label_y='<Si>(y1)')
    plot_density_curve(r=r_y2,density=magz_y2,label_y='<Si>(y2)')
    plot_density_curve(r=r_y3,density=magz_y3,label_y='<Si>(y3)')
    plot_density_curve(r=r_y4,density=magz_y4,label_y='<Si>(y4)')
    plot_density_curve(r=r_y5,density=magz_y5,label_y='<Si>(y5)')
    #
    plot_density_curve(r=r_mean,density=magz_mean,label_y='<Si>(mean)')
