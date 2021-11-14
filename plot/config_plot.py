import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import Normalize


def ax_config_plot_xyz(ax, filename, Color, LineWidth, pov="xy",rotxyz=None, xshift=0, yshift=0, zslice=None, mesh=1, bead=0, rod=0, pwlim=0, d=1):
    print(pov+" plotting",filename)
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x,y,z,ux,uy,uz,nx,ny,nz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:17]
    ns = np.transpose(data[17:])
    ens = np.array([en0, en1])
    #x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    if(pov=="zx"):
        x,y,z=z,x,y
        ux,uy,uz=uz,ux,uy

    take = z>-1000
    if zslice:
        # plot only cross section based on slice of z
        take= (z>zslice[0]) & (z<zslice[1])
    # recenter based on the sliced part
    x,y,z=x-np.average(x[take]),y-np.average(y[take]),z-np.average(z[take])

    if(rotxyz):
        rs=[x,y,z]
        us=[ux,uy,uz]
        for i in range(3):
            cosphi,sinphi = np.cos(rotxyz[i]),np.sin(rotxyz[i])
            # position
            ra=rs[(i+1)%3]*cosphi-rs[(i+2)%3]*sinphi
            rb=rs[(i+1)%3]*sinphi+rs[(i+2)%3]*cosphi
            rs[(i+1)%3],rs[(i+2)%3]=ra,rb
            # director
            ua=us[(i+1)%3]*cosphi-us[(i+2)%3]*sinphi
            ub=us[(i+1)%3]*sinphi+us[(i+2)%3]*cosphi
            us[(i+1)%3],us[(i+2)%3]=ua,ub
        x,y,z=rs
        ux,uy,uz=us

    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    alpha_xy = 0.8*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    x = x + xshift
    y = y + yshift

    # track bead index for debug
    if(0):
        for i in range(len(x)):
            pass
            ax_xy.annotate(i, (x[i], y[i]), fontsize=5)

    # beads
    if(bead):
        r = 0.5 * np.ones(len(x))
        for j in range(len(x)):
            if(take[j]):
                xc, yc, rc = x[j], y[j], r[j]
                ax.add_artist(Circle(xy=(xc, yc), linewidth=0,radius=rc,edgecolor="None",facecolor=Color, alpha=alpha_xy[j]))

    if(mesh):
        bonds = []
        for i in range(len(ns)):
            for j in range(len(ns[0])):
                if ns[i, j] != -1:
                    if(i<ns[i,j]):
                        bonds.append((i,int(ns[i, j])))
                    else:
                        bonds.append((int(ns[i, j]),i))
        bonds = set(bonds)
        for bond in bonds:
            pass
            a,b = bond
            if(take[a] or take[b]):
                ax.plot([x[a],x[b]], [y[a],y[b]], color=Color,lw=LineWidth/4,alpha=alpha_xy[a])

    # edge bond
    ecolors = ["blue","green","crimson","indigo","cyan"]
    ecolors=  ["k","k","k"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                if(take[j] or take[int(ens[i, j])]):
                    ax.plot([x[j], x[int(ens[i, j])]], [y[j], y[int(ens[i, j])]], "-",lw=LineWidth*1, color=ecolors[int(enum[j])], alpha=alpha_xy[j])

    # director
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    if(rod):
        for i in range(len(x)):
            if(take[i]):
                ax.plot([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],"-",linewidth=LineWidth,color=cmap(norm(deg[i])),solid_capstyle="round")

    if(pwlim!=0):
        deg_slct=deg>pwlim
        x_pw,y_pw,z_pw = x[deg_slct],y[deg_slct],z[deg_slct]
        ux_pw,uy_pw,uz_pw = ux[deg_slct],uy[deg_slct],uz[deg_slct]
        deg_pw = deg[deg_slct]
        for i in range(len(ux_pw)):
            ax.plot([x_pw[i]-0.5*d*sx_pw[i],x_pw[i]+0.5*d*sx_pw[i]],[y_pw[i]-0.5*d*sy_pw[i],y_pw[i]+0.5*d*sy_pw[i]],"-",linewidth=LineWidth,color=cmap(norm(deg_pw[i])),solid_capstyle="round")

    ax.set_aspect("equal")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)




def catenoid_trinoid_demo_data_get():
    foldername = "../data/Ne2/Oct23_2021"
    lf=20.0
    Cn = 5.0
    Kds = [4.2,1.8]
    Kdp = []
    fnames,povs,rotxyzs, xyshift, zslice = [],[],[],[],[]
    for i in range(len(Kds)):
        for j in range(3):
            fnames.append( foldername + "/State_N300_imod3_Ne2_lf%.1f_kar40_C00.0_karg0.0_lam5.0_Kd%.1f_q0.0_Cn%.1f.csv" % (lf,Kds[i],Cn))
            Kdp.append(Kds[i])
            povs.append("xy")
            rotxyzs.append(None)
            xyshift.append((10*j,10*i))
        # 2 crossection each
        zslice+=[(7,12),(-2.5,2.5),None]
        xyshift[3*i+2]=(25,10*i)
        povs[3*i+2]="zx"
        rotxyzs[3*i+2]=[-np.pi*3/8,np.pi*3/8,0]
    #fnames.append(fnames[-1])
    #povs.append("zx")
    #Kdp.append(Kds[-1])
    #rotxyzs.append([-np.pi*3/8,np.pi*3/8,0])
    #xyshift.append((10*2,0))
    #zslice.append(None)
    return [Kdp,fnames,povs,rotxyzs, xyshift, zslice]

def catenoid_trinoid_demo_plot(LineWidth, FontSize, LabelSize):
    print("(demo就交给我吧)")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((1, 1), (0, 0))


    folder = "../data/Ne2/Oct26_2021"
    qs = [0.0,0.8,1.6]
    for i in range(len(qs)):
        fname = folder+"/State_N300_imod3_Ne2_lf0.0_kar40_C00.0_karg0.0_lam5.0_Kd6.0_q%.1f_Cn6.0.csv"%qs[i]
        ax_config_plot_xyz(axcfg, fname, "gray", LineWidth, pov="zx", rotxyz=(0,np.pi/8,0),xshift=20*i,yshift=0, mesh=1, bead=0,rod=0,d=0.7)
        ax_config_plot_xyz(axcfg, fname, "gray", LineWidth, pov="zx", rotxyz=(0,np.pi/8,0),xshift=20*i,yshift=-20, mesh=1, bead=0,rod=1,d=0.7)
        axcfg.text(20*i,-33,r"$%.1f$"%qs[i],fontsize=FontSize)
    axcfg.text(-20,-33,r"$k_c\sqrt{\frac{\epsilon_{LL}}{C}}=$",fontsize=FontSize)
    i=3
    qs.append(2.0)
    #folder = "../data/Ne2/Oct23_2021"
    fname = folder+"/State_N300_imod3_Ne3_lf0.0_kar40_C00.0_karg0.0_lam5.0_Kd6.0_q2.0_Cn6.0.csv"
    ax_config_plot_xyz(axcfg, fname, "gray", LineWidth, pov="xy", rotxyz=(0,0,0),xshift=20*i,yshift=0, mesh=1, bead=0,rod=0,d=0.7)
    ax_config_plot_xyz(axcfg, fname, "gray", LineWidth, pov="xy", rotxyz=(0,0,0),xshift=20*i,yshift=-20, mesh=1, bead=0,rod=1,d=0.7)
    axcfg.text(20*i,-33,r"$%.1f$,(3 edges)"%(qs[i]*np.sqrt(6/6)),fontsize=FontSize)

    axcfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)

    plt.tight_layout(pad=0.1)
    plt.savefig("catenoid_demo.pdf",format="pdf")
