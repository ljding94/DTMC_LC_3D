import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.special import expit

def ax_config_plot_xyz(ax, filename, Color, LineWidth, pov="xy",rotxyz=None, xshift=0, yshift=0, zslice=None, hide_behind=None, mesh=1, twistcolor=0 ,bead=0, rod=0, pwlim=0, d=1, scale=1):
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

    x,y,z = scale*x,scale*y,scale*z

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

    if(hide_behind):
        z_cut = z_min + (z_max-z_min)*hide_behind
        take = take & (z>z_cut)


    #alpha_xy = np.maximum(1.0*(z-z_min+0.1)/(z_max-z_min+0.1)-0.2,0)
    alpha_xy = 0.8*(z-z_min+0.01)/(z_max-z_min+0.01)+0.1
    #alpha_rod = 0.8*(z-z_min+0.01)/(z_max-z_min+0.01)+0.1
    al = 1.2
    #alpha_xy = expit(1.0*(z-0.5*z_min))
    alpha_rod = expit(al*(z-0.5*z_min))
    #alpha_rod = np.maximum(1.0*(z-z_min+0.01)/(z_max-z_min+0.01)-0.1,0)
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
                ax.add_artist(Circle(xy=(xc, yc), linewidth=0,radius=rc,edgecolor="None",facecolor=Color, alpha=0.8*alpha_xy[j]))

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
        if(twistcolor):
            normt=Normalize(vmin=0,vmax=0.25*np.pi)
            cmapt = cm.get_cmap("ocean")
        for bond in bonds:
            a,b = bond
            if(take[a] and take[b]):
                if(twistcolor):
                    twist = twist_cal([x[a],y[a],z[a]],[x[b],y[b],z[b]],[ux[a],uy[a],uz[a]],[ux[b],uy[b],uz[b]])
                    ax.plot([x[a],x[b]], [y[a],y[b]], color=cmapt(normt(twist)),lw=LineWidth/4,alpha=alpha_xy[a])
                else:
                    ax.plot([x[a],x[b]], [y[a],y[b]], color=Color,lw=LineWidth/4,alpha=alpha_xy[a])


    # director
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    if(rod):
        for i in range(len(x)):
            if(take[i]):
                ax.plot([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],"-",linewidth=LineWidth,color=cmap(norm(deg[i])),solid_capstyle="round",alpha=alpha_rod[i])

    if(pwlim!=0):
        deg_slct=deg>pwlim
        #deg_slct=np.sqrt(un2)<pwlim
        x_pw,y_pw,z_pw = x[deg_slct],y[deg_slct],z[deg_slct]
        ux_pw,uy_pw,uz_pw = ux[deg_slct],uy[deg_slct],uz[deg_slct]
        deg_pw = deg[deg_slct]
        alpha_xy_pw=alpha_xy[deg_slct]
        r = 0.5 * np.ones(len(x_pw))
        for i in range(len(ux_pw)):
            if(take[i]):
                xc, yc, rc = x_pw[i], y_pw[i], r[i]
                ax.add_artist(Circle(xy=(xc, yc), linewidth=0,radius=rc,edgecolor="None",facecolor=Color, alpha=0.8*alpha_rod[i]))
                ax.plot([x_pw[i]-0.5*d*ux_pw[i],x_pw[i]+0.5*d*ux_pw[i]],[y_pw[i]-0.5*d*uy_pw[i],y_pw[i]+0.5*d*uy_pw[i]],"-",linewidth=LineWidth,color=cmap(norm(deg_pw[i])),solid_capstyle="round",alpha=alpha_rod[i])

    # edge bond, plot edge lastly, so it's always shown
    ecolors = ["blue","green","crimson","indigo","cyan"]
    ecolors=  ["k","k","k"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                if(take[j] or take[int(ens[i, j])]):
                    ax.plot([x[j], x[int(ens[i, j])]], [y[j], y[int(ens[i, j])]], "-",lw=LineWidth*1, color=ecolors[int(enum[j])], alpha=alpha_xy[j])

    ax.set_aspect("equal")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)


def twist_cal(r1,r2,u1,u2):
    r12 = np.array(r2)-np.array(r1)
    r12/=np.sqrt(np.dot(r12,r12))
    ans = np.dot(np.cross(u1,u2),r12)*np.dot(u1,u2)
    return ans

def config_demo(LineWidth, FontSize, LabelSize):
    print("👌")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.4))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax = plt.subplot2grid((1, 1), (0, 0))
    msize=4
    fname = "../data/Ne2/Feb_2022/Feb28_2022/State_N300_imod3_Ne2_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_q1.4_Cn4.0.csv"
    rot = (1.6,0.0,0)
    # configuration plot
    ax_config_plot_xyz(ax, fname, "gray", LineWidth,rotxyz=(1.0,0.0,0), xshift=0, mesh=1, bead=0,rod=0, d=0.8)
    x0=4
    ax.text(0-x0, -10, r"(a)", fontsize=FontSize)
    ax_config_plot_xyz(ax, fname, "gray", LineWidth, rotxyz=(1.4,0.0,0),xshift=17, mesh=1, bead=1,rod=1, d=0.8)
    ax.text(17-x0, -10, r"(b)", fontsize=FontSize)
    ax_config_plot_xyz(ax, fname, "gray", LineWidth, rotxyz=(1.8,0.0,0),xshift=34, mesh=0, bead=0, rod=0, d=0.8,pwlim=np.pi/3)
    ax.text(34-x0, -10, r"(c)", fontsize=FontSize)


    #cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),ax=ax,ticks=[0,0.25*np.pi,0.5*np.pi])
    #cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"],fontsize=FontSize)
    cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),ax=ax,ticks=[0,np.pi/6,np.pi/3,np.pi/2])
    cbar.ax.set_yticklabels([r"$0$",r"$\pi/6$",r"$\pi/3$",r"$\pi/2$"],fontsize=FontSize)
    cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    cbar.ax.set_title(r"$\arccos{|\vu{u}\cdot\vu{n}|}$",fontsize=FontSize)

    ax.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    #x1, y1 = 0.9, 1.06
    #axcfg.text(x1,y1, r"(a)", fontsize=FontSize,transform=axuz.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/config_demo.pdf",format="pdf")


def init_config_demo(LineWidth, FontSize, LabelSize):
    print("plotting figures for initial shapes")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.4))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax = plt.subplot2grid((1, 1), (0, 0))
    msize=4
    rot = (1.6,0.0,0)
    # configuration plot
    ax_config_plot_xyz(ax, "init_config/flat.csv", "gray", LineWidth,rotxyz=(0,0.0,0), xshift=0, mesh=1, bead=0,rod=0,scale=0.5)
    x0=4
    ax.text(0-x0, -10, r"(a)", fontsize=FontSize)
    ax_config_plot_xyz(ax, "init_config/cylinder.csv", "gray", LineWidth, pov="zx", rotxyz=(0.2,0.4,0),xshift=17, mesh=1, bead=0,rod=0)
    ax.text(17-x0, -10, r"(b)", fontsize=FontSize)
    ax_config_plot_xyz(ax, "init_config/vesicle.csv", "gray", LineWidth, rotxyz=(0.0,0.6,0),xshift=34, mesh=1, bead=0, rod=0)
    ax.text(34-x0, -10, r"(c)", fontsize=FontSize)

    ax.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    #x1, y1 = 0.9, 1.06
    #axcfg.text(x1,y1, r"(a)", fontsize=FontSize,transform=axuz.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/init_config_demo.pdf",format="pdf")

def cylinder_config_demo(LineWidth, FontSize, LabelSize):
    print("plotting figures for cylinder config demo")

    lf = 25
    foldername = "../data/Ne2/May12_2022"
    fname = foldername + "/State_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.0_Cn4.0.csv"%lf

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax = plt.subplot2grid((1, 1), (0, 0))
    msize=4
    dely=8
    ax_config_plot_xyz(ax, fname, "gray", LineWidth, pov="zx", mesh=1, bead=0,rod=0, d=0.8)

    ax.text(lf/2+1,-2, r"(a)",fontsize=FontSize)
    ax_config_plot_xyz(ax, fname, "gray", LineWidth, pov="zx", yshift=-dely, mesh=1, bead=1,rod=1, d=0.8)
    ax.text(lf/2+1,-2-dely, r"(b)",fontsize=FontSize)

    delb=5
    ax.plot([-lf/2,-lf/2],[-dely-delb,delb], color="black", ls="--", lw = LineWidth)
    ax.plot([lf/2,lf/2],[-dely-delb,delb], color="black", ls="--", lw = LineWidth)
    ax.text(-lf/2,-dely-delb-2, r"$z=0$",fontsize=FontSize)
    ax.text(lf/2,-dely-delb-2, r"$z=l_f$",fontsize=FontSize)



    plt.tight_layout(pad=0.1)
    #plt.show()
    plt.savefig("figures/elongation_config_demo.pdf",format="pdf")




def twist_config_demo():
    pass
    # configuration comparison, director color map, and twist color on bond



# not paper related

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

def catenoid_demo_3d_plot():
    #TODO: modify the following for the catenoid config accordingly
    print("plotting configuration demo plot for slide")
    filename="../data/Ne2/Feb2_2022/State_N300_imod3_Ne3_lf0.0_kar100_C00.0_karg0.0_lam5.0_Kd5.0_q2.0_Cn5.0.csv"
    mesh,rod=1,0
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    #x,y,z,sx,sy,sz,enum, en0, en1 = data[5:14]

    x,y,z,sx,sy,sz,nx,ny,nz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:17]
    #x,y,z,sx,sy,sz,phi,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:15]
    ns = np.transpose(data[17:])
    # just for illustrating surface normal
    #x,y,z,sx,sy,sz,nx,ny,nz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:17]
    #ns = np.transpose(data[17:])

    #x,y,z,sx,sy,sz,dA,d2H,ds,dAK,un2,is_cnp,enum, en0, en1 = data[:15]
    d=1.5
    #sx,sy,sz=d*sx,d*sy,d*sz
    x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    x_min, x_max = np.min(x),np.max(x)
    y_min, y_max = np.min(y),np.max(y)
    z_min, z_max = np.min(z),np.max(z)
    max_range_half = max([x_max-x_min,y_max-y_min,z_max-z_min])*0.5
    alpha_xy = 0.9*(z-z_min+0.1)/(z_max-z_min+0.1)+0.1
    alpha_zx = 0.9*(y-y_min+0.1)/(y_max-y_min+0.1)+0.1

    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(5, 5))
    ax = plt.axes(projection="3d")

    #for i in range(len(x)):
    #    ax.text(x[i],y[i],z[i],i)

    if(mesh):
        for i in range(len(ns)):
            for j in range(len(ns[0])):
                if ns[i, j] != -1:
                    pass
                    ax.plot3D([x[i], x[int(ns[i, j])]], [y[i], y[int(ns[i, j])]], [z[i], z[int(ns[i, j])]], "-",  color="silver")

    ecolors = ["blue","green","crimson","indigo","cyan","black"]
    #ecolors = ["blue","purple","yellow","cyan","red","green","black"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax.plot3D([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], [
                    z[j], z[int(ens[i, j])]], "-",linewidth=2, color=ecolors[int(enum[j])], alpha=0.5)
    cmap = cm.get_cmap("jet_r")

    # director
    x_ave,y_ave,z_ave = np.average(x),np.average(y),np.average(z)
    #deg = np.arccos(np.absolute(snu))
    #deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    deg = np.arccos(np.sqrt(un2))
    if(rod):
        for i in range(len(sx)):
            ax.plot3D([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],[z[i]-0.5*d*sz[i],z[i]+0.5*d*sz[i]],"-",color=cmap(norm(deg[i])))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        #cbar=plt.colorbar(sm, ticks=[0,0.25*np.pi,0.5*np.pi])
        #cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"])
        #cbar.ax.set_title(r"$\arccos|u\cdot n|$")

    ax.set_xlim(-max_range_half,max_range_half)
    ax.set_ylim(-max_range_half,max_range_half)
    ax.set_zlim(-max_range_half,max_range_half)
    ax.set_axis_off()

    for ii in range(20):
        ax.view_init(elev=30., azim=36/2*ii)
        plt.tight_layout()
        #plt.savefig("figures/config3d_demo/demo_rod_ani%d"%ii+".png",dpi=300)
        plt.savefig("figures/catenoid3d_demo/demo_ani%d"%ii+".png",dpi=300)
    plt.close()