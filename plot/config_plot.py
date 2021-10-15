import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import Normalize


def ax_config_plot_xy(ax, filename, Color,Line2D, rotation=None, xshift=0, yshift=0, mesh=1, bead=0, rod=0, pwlim=0, d=1):
    print("ax plotting",filename)
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x,y,z,ux,uy,uz,nx,ny,nz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:17]
    ns = np.transpose(data[17:])
    ens = np.array([en0, en1])
    x,y,z=x-np.average(x),y-np.average(y),z-np.average(z)
    if(rotation):
        # rotate beads
        phi,theta=rotation
        xr=x*np.cos(phi)-y*np.sin(phi)
        yr = np.cos(theta)*(x*np.sin(phi)+y*np.cos(phi))-z*np.sin(theta)
        zr = np.sin(theta)*(x*np.sin(phi)+y*np.cos(phi))+z*np.cos(theta)
        x,y,z=xr,yr,zr
        # rotate directors
        uxr=ux*np.cos(phi)-uy*np.sin(phi)
        uyr = np.cos(theta)*(ux*np.sin(phi)+uy*np.cos(phi))-uz*np.sin(theta)
        uzr = np.sin(theta)*(ux*np.sin(phi)+uy*np.cos(phi))+uz*np.cos(theta)
        ux,uy,uz=uxr,uyr,uzr

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
            ax.plot([x[a],x[b]], [y[a],y[b]], color=Color,lw=LineWidth/4,alpha=alpha_xy[a])

    # edge bond
    ecolors = ["blue","green","crimson","indigo","cyan"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax.plot([x[j], x[int(ens[i, j])]], [
                    y[j], y[int(ens[i, j])]], "-",lw=LineWidth/2, color=ecolors[int(enum[j])], alpha=alpha_xy[j])

    # director
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    if(rod):
        ftail+="_rod"
        for i in range(len(x)):
            ax_xy.plot([x[i]-0.5*d*sx[i],x[i]+0.5*d*sx[i]],[y[i]-0.5*d*sy[i],y[i]+0.5*d*sy[i]],"-",linewidth=LineWidth,color=cmap(norm(deg[i])),solid_capstyle="round")

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
