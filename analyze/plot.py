import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import Normalize
from analyze import *
from scipy import optimize


def find_nu(ux, uy, uz):
    Q = np.array([[1.5 * np.average(ux * ux) - 0.5, np.average(ux * uy), 1.5 * np.average(ux * uz)], [0, 1.5 * np.average(uy * uy) - 0.5, np.average(uy * uz)], [0, 0, 1.5 * np.average(uz * uz) - 0.5]])
    Q[1, 0] = Q[0, 1]
    Q[2, 0] = Q[0, 2]
    Q[2, 1] = Q[1, 2]
    w, v = np.linalg.eig(Q)
    w_max = np.max(w)
    for i in range(len(w)):
        if w[i] == w_max:
            return np.transpose(v)[i]


def mean_filter(heat, ns):
    # smooth the heat map by averaging with neighbors bases on ns
    mean_heat = []
    for i in range(len(heat)):
        mean_heat.append(heat[i])
        nei_count = 0
        for j in range(len(ns[0])):
            if ns[i, j] != -1:
                mean_heat[i] += heat[int(ns[i, j])]
                nei_count += 1
        mean_heat[i] /= 1 + nei_count
    return mean_heat


def config_plot_xyz(filename, mesh=0, rod=1, piwall=0, phicolor=0, cvt_map="", cmap_smooth=0, tag="", Format="pdf", lim=15, fix_index=None):
    print("plotting", filename)
    ftail = "_xyz"
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)

    x, y, z, sx, sy, sz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    # x,y,z,sx,sy,sz,nx,ny,nz,phi,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:18]
    ns = np.transpose(data[17:])
    # ns = np.transpose(data[18:])
    # sx,sy,sz=d*sx,d*sy,d*sz
    # x,y,z,sx,sy,sz, enum, en0, en1 = data[5:14]
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    z_min, z_max = np.min(z), np.max(z)
    alpha_xy = 0.8 * (z - z_min + 0.1) / (z_max - z_min + 0.1) + 0.1
    alpha_zx = 0.8 * (y - y_min + 0.1) / (y_max - y_min + 0.1) + 0.1
    # ns = np.transpose(data[14:])
    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(10, 5))
    # ax_xy = fig.add_subplot(111, aspect="equal")
    ax_xy = fig.add_subplot(121, aspect="equal")
    ax_zx = fig.add_subplot(122, aspect="equal")
    # beads
    # phi field
    if phicolor:
        phicmap = cm.get_cmap("PuOr")
        phinorm = Normalize(vmin=-1, vmax=1)
        ax_xy.scatter(x, y, c=phi, marker="o", cmap=phicmap, norm=phinorm)
        # ax_xy.scatter(x[phi==1],y[phi==1],marker="o",facecolor="None",edgecolor="black")
        ax_zx.scatter(z, x, c=phi, marker="o", cmap=phicmap, norm=phinorm)
        # ax_zx.scatter(z[phi==1],x[phi==1],marker="o",facecolor="None",edgecolor="black")
        phism = plt.cm.ScalarMappable(cmap=phicmap, norm=phinorm)
        phism.set_array([])
        phicbar = plt.colorbar(phism, ticks=[-1, -0.5, 0, 0.5, 1])
        phicbar.ax.set_title(r"$\phi$")
    # bulk bond

    # track bead ind #
    if 0:
        for i in range(len(x)):
            pass
            ax_xy.annotate(i, (x[i], y[i]), fontsize=5)

    if mesh:
        bonds = []
        for i in range(len(ns)):
            for j in range(len(ns[0])):
                if ns[i, j] != -1:
                    if i < ns[i, j]:
                        bonds.append((i, int(ns[i, j])))
                    else:
                        bonds.append((int(ns[i, j]), i))
        bonds = set(bonds)
        for bond in bonds:
            pass
            a, b = bond
            ax_xy.plot([x[a], x[b]], [y[a], y[b]], color="silver", alpha=alpha_xy[a])
            ax_zx.plot([z[a], z[b]], [x[a], x[b]], color="silver", alpha=alpha_xy[a])
        # color mapping
    # edge bond
    ecolors = ["blue", "green", "crimson", "indigo", "cyan"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax_xy.plot([x[j], x[int(ens[i, j])]], [y[j], y[int(ens[i, j])]], "-", linewidth=0.5, color=ecolors[int(enum[j])], alpha=alpha_xy[j])
                ax_zx.plot([z[j], z[int(ens[i, j])]], [x[j], x[int(ens[i, j])]], "-", linewidth=2, color=ecolors[int(enum[j])], alpha=alpha_zx[j])
    # spin vector
    """'
    nu=find_nu(sx,sy,sz)
    #nu=[0,0,1]
    x_ave,y_ave,z_ave = np.average(x),np.average(y),np.average(z)
    D_ave = 3
    ax_xy.plot([x_ave-D_ave*nu[0],x_ave+D_ave*nu[0]],[y_ave-D_ave*nu[1],y_ave+D_ave*nu[1]],"-",linewidth=3.0,color="k")
    ax_zx.plot([z_ave-D_ave*nu[2],z_ave+D_ave*nu[2]],[x_ave-D_ave*nu[0],x_ave+D_ave*nu[0]],"-",linewidth=3.0,color="k")
    snu=sx*nu[0]+sy*nu[1]+sz*nu[2]
    snu=sz
    """
    # print(np.sqrt(un2)-np.absolute(snu))
    # deg = np.arccos(np.absolute(snu))
    norm = Normalize(vmin=0, vmax=0.5 * np.pi)
    deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    if rod:
        ftail += "_rod"
        for i in range(len(x)):
            d = 2
            line = "-"
            ax_xy.plot([x[i] - 0.5 * d * sx[i], x[i] + 0.5 * d * sx[i]], [y[i] - 0.5 * d * sy[i], y[i] + 0.5 * d * sy[i]], linestyle=line, linewidth=1.5, color=cmap(norm(deg[i])))
            ax_zx.plot([z[i] - 0.5 * d * sz[i], z[i] + 0.5 * d * sz[i]], [x[i] - 0.5 * d * sx[i], x[i] + 0.5 * d * sx[i]], linestyle=line, linewidth=1.5, color=cmap(norm(deg[i])))
        # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        # sm.set_array([])
        # cbar=plt.colorbar(sm, ticks=[0,0.25*np.pi,0.5*np.pi])
        # cbar.ax.set_yticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"])
    if piwall:
        deg_slct = deg > 0.4 * np.pi
        x_pw, y_pw, z_pw = x[deg_slct], y[deg_slct], z[deg_slct]
        sx_pw, sy_pw, sz_pw = sx[deg_slct], sy[deg_slct], sz[deg_slct]
        deg_pw = deg[deg_slct]
        for i in range(len(sx_pw)):
            d = 2
            line = "-"
            ax_xy.plot([x_pw[i] - 0.5 * d * sx_pw[i], x_pw[i] + 0.5 * d * sx_pw[i]], [y_pw[i] - 0.5 * d * sy_pw[i], y_pw[i] + 0.5 * d * sy_pw[i]], linestyle=line, linewidth=1.5, color=cmap(norm(deg_pw[i])))
            ax_zx.plot([z_pw[i] - 0.5 * d * sz_pw[i], z_pw[i] + 0.5 * d * sz_pw[i]], [x_pw[i] - 0.5 * d * sx_pw[i], x_pw[i] + 0.5 * d * sx_pw[i]], linestyle=line, linewidth=1.5, color=cmap(norm(deg_pw[i])))

    # plot fixed bead (if have)
    if fix_index:
        ax_xy.plot([x[fix_index[0]], x[fix_index[1]]], [y[fix_index[0]], y[fix_index[1]]], marker="o", linestyle="None", color="purple")
        ax_zx.plot([z[fix_index[0]], z[fix_index[1]]], [x[fix_index[0]], x[fix_index[1]]], marker="o", linestyle="None", color="purple")
    ax_xy.set_xlim(x_min - 2, x_max + 2)
    ax_xy.set_ylim(y_min - 2, y_max + 2)
    # ax_xy.set_title("XY  "+tag+" smooth=%d"%cmap_smooth, fontsize=25)
    ax_zx.set_xlim(z_min - 2, z_max + 2)
    ax_zx.set_ylim(x_min - 2, x_max + 2)
    ax_zx.set_title("ZX")
    ax_xy.legend(title=tag)
    # plt.savefig(filename[:-4] + "_xy.pdf", dpi=300, format="pdf")
    plt.savefig(filename[:-4] + ftail + "." + Format, dpi=100, format=Format, bbox_inches="tight", transparent=False)
    plt.close()


def config_plot_xyz_seq(filename, Seq):
    for i in range(Seq):
        config_plot_xyz(filename[:-4] + "_%d.txt" % i, Format="png")


def config_plot3D(filename, mesh=0, rod=0, piwall=0, phicolor=0, fnormal=0, cvt_map="", cmap_smooth=0):
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    # x,y,z,sx,sy,sz,enum, en0, en1 = data[5:14]

    x, y, z, sx, sy, sz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    # x,y,z,sx,sy,sz,phi,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:15]
    ns = np.transpose(data[17:])
    # just for illustrating surface normal
    # x,y,z,sx,sy,sz,nx,ny,nz,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:17]
    # ns = np.transpose(data[17:])

    # x,y,z,sx,sy,sz,dA,d2H,ds,dAK,un2,is_cnp,enum, en0, en1 = data[:15]
    d = 1
    # sx,sy,sz=d*sx,d*sy,d*sz
    x, y, z = x - np.average(x), y - np.average(y), z - np.average(z)
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    z_min, z_max = np.min(z), np.max(z)
    max_range_half = max([x_max - x_min, y_max - y_min, z_max - z_min]) * 0.5
    alpha_xy = 0.9 * (z - z_min + 0.1) / (z_max - z_min + 0.1) + 0.1
    alpha_zx = 0.9 * (y - y_min + 0.1) / (y_max - y_min + 0.1) + 0.1

    ens = np.array([en0, en1])
    fig = plt.figure(figsize=(5, 5))
    ax = plt.axes(projection="3d")

    # for i in range(len(x)):
    #    ax.text(x[i],y[i],z[i],i)

    if mesh:
        for i in range(len(ns)):
            for j in range(len(ns[0])):
                if ns[i, j] != -1:
                    pass
                    ax.plot3D([x[i], x[int(ns[i, j])]], [y[i], y[int(ns[i, j])]], [z[i], z[int(ns[i, j])]], "-", color="silver")

    ecolors = ["blue", "green", "crimson", "indigo", "cyan", "black"]
    # ecolors = ["blue","purple","yellow","cyan","red","green","black"]
    for i in range(len(ens)):
        for j in range(len(en0)):
            if ens[i, j] != -1:
                ax.plot3D([x[j], x[int(ens[i, j])]], [y[j], y[int(ens[i, j])]], [z[j], z[int(ens[i, j])]], "-", linewidth=2, color=ecolors[int(enum[j])], alpha=0.5)
    cmap = cm.get_cmap("jet_r")
    if cvt_map == "Mean":
        # ftail+="_mmap"
        # norm=Normalize(vmin=0,vmax=0.5*np.pi)
        # heat = dA*d2H
        heat = d2H
        for m in range(cmap_smooth):
            heat = mean_filter(heat, ns)
        ax.scatter3D(x, y, z, c=cmap(heat), s=10)
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_array([])
        cbar = plt.colorbar(sm)

    # director
    nu = find_nu(sx, sy, sz)
    print(nu, np.dot(nu, nu))
    x_ave, y_ave, z_ave = np.average(x), np.average(y), np.average(z)
    D_ave = 2
    # ax.plot3D([x_ave-D_ave*nu[0],x_ave+D_ave*nu[0]],[y_ave-D_ave*nu[1],y_ave+D_ave*nu[1]],[z_ave-D_ave*nu[2],z_ave+D_ave*nu[2]],"-",linewidth=3.0,color="k")
    snu = sx * nu[0] + sy * nu[1] + sz * nu[2]
    # deg = np.arccos(np.absolute(snu))
    # deg = np.arccos(np.sqrt(un2))
    cmap = cm.get_cmap("jet_r")
    norm = Normalize(vmin=0, vmax=0.5 * np.pi)
    deg = np.arccos(np.sqrt(un2))
    if rod:
        for i in range(len(sx)):
            ax.plot3D([x[i] - 0.5 * d * sx[i], x[i] + 0.5 * d * sx[i]], [y[i] - 0.5 * d * sy[i], y[i] + 0.5 * d * sy[i]], [z[i] - 0.5 * d * sz[i], z[i] + 0.5 * d * sz[i]], "-", color=cmap(norm(deg[i])))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ticks=[0, 0.25 * np.pi, 0.5 * np.pi])
        cbar.ax.set_yticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$"])
        cbar.ax.set_title(r"$\arccos|u\cdot n|$")
    if piwall:
        deg_slct = deg > 0.4 * np.pi
        x_pw, y_pw, z_pw = x[deg_slct], y[deg_slct], z[deg_slct]
        sx_pw, sy_pw, sz_pw = sx[deg_slct], sy[deg_slct], sz[deg_slct]
        deg_pw = deg[deg_slct]
        # print("deg_pw[0]",deg_pw[0])
        for i in range(len(sx_pw)):
            ax.plot3D([x_pw[i] - 0.5 * d * sx_pw[i], x_pw[i] + 0.5 * d * sx_pw[i]], [y_pw[i] - 0.5 * d * sy_pw[i], y_pw[i] + 0.5 * d * sy_pw[i]], [z_pw[i] - 0.5 * d * sz_pw[i], z_pw[i] + 0.5 * d * sz_pw[i]], "-", color=cmap(norm(deg_pw[i])))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ticks=[0, 0.25 * np.pi, 0.5 * np.pi])
        cbar.ax.set_yticklabels([r"$0$", r"$\pi/4$", r"$\pi/2$"])
        cbar.ax.set_title(r"$\arccos|u\cdot n|$")

    if fnormal:
        for i in range(len(sx)):
            ax.plot3D([x[i], x[i] + d * nx[i]], [y[i], y[i] + d * ny[i]], [z[i], z[i] + d * nz[i]], "k-")

    ax.set_xlim(-max_range_half, max_range_half)
    ax.set_ylim(-max_range_half, max_range_half)
    ax.set_zlim(-max_range_half, max_range_half)
    # ax.scatter3D(x, y, z, s=[1 for i in range(len(x))])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_axis_off()
    plt.show()
    # plt.savefig(filename[:-4] + "_3D.png", dpi=300)
    plt.close()


def config_nu2_dis(filename, bin_num=20):
    ftail = "_dis"
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x, y, z, sx, sy, sz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:14]
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.hist(un2, bins=bin_num, histtype="step", density=True)
    ax1.hist(dA * d2H * d2H, bins=bin_num, histtype="step", density=True)
    ax2.scatter(un2, dA * d2H * d2H)
    plt.show()


def autocorrelation_plot(rho, tau_int, savefile):
    t = np.linspace(0, 1000, 1000)
    plt.figure()
    plt.plot(range(1000), rho[:1000], "d")
    plt.plot(t, np.exp(-t / tau_int), "--")
    plt.savefig(savefile, format="pdf", transparent=True)
    plt.close()


def O_MCstep_plot(filename, thermN, Ne):
    data = np.loadtxt(filename, skiprows=14, delimiter=",", unpack=True)
    E = data[0]
    Les = data[1 : 1 + Ne]
    Lasym = eff_Le_sym_cal(Les)
    (
        IdA,
        I2H,
        I2H2,
        I2H2dis,
        IK,
        Tp2uu,
        Tuuc,
        Bond_num,
        Tun2,
        Tuz2,
        Tlb,
    ) = data[1 + Ne : 12 + Ne]
    # Eu  = data[12+Ne]
    p2uu = Tp2uu / Bond_num
    uuc = Tuuc / Bond_num
    lb = Tlb / Bond_num
    # Euave = np.average(Eu)
    # wnu = np.exp(Eu-Euave)  # biased we

    ppi = 72
    plt.figure()
    plt.rc("text")
    na = 5
    fig, axs = plt.subplots(na, 1, figsize=(246 / ppi * 1, 246 / ppi * 0.5 * na), sharex=True)  # , sharex=True
    mcs = np.arange(1, len(E) + 0.1, 1)
    axs[0].plot(mcs, E, "-")
    axs[0].set_ylabel("E")
    axs[1].plot(mcs[thermN:], E[thermN:], "-")
    axs[1].set_ylabel("E")
    for i in range(len(Les)):
        axs[2].plot(mcs, Les[i], "-")
        print(Les[i])
    axs[2].set_ylabel("Le")
    axs[3].plot(mcs, Lasym, "-")
    axs[3].set_ylabel("Lasym")
    axs[4].plot(mcs, I2H2, "-")
    axs[4].set_ylabel("I2H2")
    axs[-1].set_xlabel("MC steps")
    plt.tight_layout()
    plt.savefig(filename[:-4] + "_MCstep.png")
    plt.close()


def uz_distribution_plot(filename):
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    # x,y,z,sx,sy,sz,enum, en0, en1 = data[5:14]
    x, y, z, sx, sy, sz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    # x,y,z,sx,sy,sz,phi,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:15]
    ns = np.transpose(data[17:])
    plt.figure()
    plt.hist(np.arccos(sz), bins=50, density=1)
    plt.xlabel("arccos(uz)")
    plt.ylabel("density")
    plt.show()
    plt.close()


# TODO: finished alpha distribution function


def tan_alpha_distribution_plot(filename):
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    # x,y,z,sx,sy,sz,enum, en0, en1 = data[5:14]
    x, y, z, ux, uy, uz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    ns = np.transpose(data[17:])
    # x,y,z,sx,sy,sz,phi,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:15]
    # find center (x0,y0) for the cylinderical part
    x, y, z = x - np.average(x), y - np.average(y), z - np.average(z)
    x_uni, y_uni = x / (np.sqrt(x ** 2 + y ** 2)), y / (np.sqrt(x ** 2 + y ** 2))

    tan_alpha = uz / (-y_uni * ux + x_uni * uy)  # u*z/u*t_phi #t_phi = (-y,x)

    fig, axs = plt.subplots(2, 1, figsize=(4, 8))
    axs[0].scatter(x, y)
    axs[1].hist(np.arctan(tan_alpha), bins=50, density=1)
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    axs[1].set_xlabel(r"$\alpha(\gamma=1)$")
    # axs[1].set_xlim(-5,5)
    axs[1].set_ylabel("density")
    plt.show()
    plt.close()





# pi wall wrapping angle analysis

def test_func_sin(phi,a,b,c):
    return a+b*np.sin(2*(phi-c))

def test_func_sin1(phi,c):
    return 0.5+0.5*np.sin(2*(phi-c))

def test_func_abx(x,a,b):
    return a+b*np.array(x)

def un2_fit_phi(phi,un2,phi0_bound):
    params, params_covariance = optimize.curve_fit(test_func_sin, phi, un2, p0=[0.5, 0.5, np.average(phi0_bound)],bounds=((0.45,0.45,phi0_bound[0]),(0.55,0.55,phi0_bound[1])))
    #params, params_covariance = optimize.curve_fit(test_func_sin1, phi, un2, p0=phi0_bound[0],bounds=(phi0_bound[0],phi0_bound[1]))
    print("params: ", params)
    #print("params_covariance: ",params_covariance)
    return params

def tilt_slice_distri_plot(filename):
    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x, y, z, ux, uy, uz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    ns = np.transpose(data[17:])
    # x,y,z,sx,sy,sz,phi,dA,d2H,ds,dAK,un2,enum, en0, en1 = data[:15]
    # find center (x0,y0) for the cylinderical part
    x, y, z = x - np.average(x), y - np.average(y), z - np.average(z)
    r = np.sqrt(x**2+y**2)
    zmin, zmax = np.min(z), np.max(z)
    zmin, zmax = zmin + (zmax - zmin) * 1 / 4, zmin + (zmax - zmin) * 3 / 4
    phi = np.arctan2(y, x)

    fig = plt.figure(figsize=(9, 6))
    ax0 = plt.subplot2grid((2, 3), (0, 0))
    ax1 = plt.subplot2grid((2, 3), (1, 0))
    ax2 = plt.subplot2grid((2, 3), (0, 1), rowspan=2)
    ax3 = plt.subplot2grid((2, 3), (0, 2), rowspan=2)
    axs = [ax0,ax1,ax2,ax3]

    nbin = 5
    colors = ["red","blue","green","tomato","black","purple"]
    z_mean = []
    R_mean = []
    phi0 = []
    for i in range(nbin):
        zl, zr = zmin + (zmax - zmin) * i / nbin, zmin + (zmax - zmin) * (i + 1) / nbin
        select = np.logical_and(z > zl, z <= zr)
        z_mean.append(np.average(z[select]))
        R_mean.append(np.average(r[select]))
        axs[0].scatter(x[select], y[select],color=colors[i])
        axs[1].hist(un2[select], histtype="step", label=r"$z\in[%.1f,%.1f)$" % (zl, zr), bins=15, color = colors[i])

        if(i==0):
            phi0_bound = [-np.pi/2,np.pi/2]
        if (1):
            un2_sort = un2[select][phi[select].argsort()]
            phi_sort = np.sort(phi[select])
            #print("phi_sort", phi_sort)
            #print("un2_sort", un2_sort)
            # axs[2].scatter(phi[select],un2[select])
            axs[2].plot(phi_sort, un2_sort+i, "o:",markersize=5, alpha=0.7, color = colors[i])
            # plot fit sin curve
            print("phi0_bound",phi0_bound)
            para_fit = un2_fit_phi(phi_sort,un2_sort,phi0_bound)
            phi0_bound = [para_fit[2]-1,para_fit[2]+1]

            axs[2].plot(phi_sort,test_func_sin(phi_sort,para_fit[0],para_fit[1],para_fit[2])+i, color = colors[i], label = r"$\phi_0=%.1f$"%para_fit[2])
            #axs[2].plot(phi_sort,test_func_sin1(phi_sort,para_fit[0])+i*0.5, color = colors[i], label = r"$\phi_0=%.1f$"%para_fit[0])
            phi0.append(para_fit[2])
            #print("para_fit ", para_fit)
    z_r_mean = np.array(z_mean)/np.array(R_mean)
    axs[3].plot(z_r_mean, phi0, "o", markersize=8, mfc="None")
    params, pcov = optimize.curve_fit(test_func_abx,z_r_mean, phi0, p0=[0, 0.1])
    perr = np.sqrt(np.diag(pcov))
    axs[3].plot(z_mean,test_func_abx(z_r_mean,params[0],params[1]),"-",label=r"$\phi_0 = %.1f+%.3f z/r$"%(params[0],params[1]))
    print(filename, params[1])
    # axs[2].scatter(phi,x)

    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    axs[1].set_xlabel(r"$(u\cdot n)^2$")
    # axs[1].set_xlim(-5,5)
    axs[1].set_ylabel("density")
    axs[1].legend()

    axs[2].set_xlabel(r"$\phi = \arctan(y/x)$")
    axs[2].set_ylabel(r"$(u\cdot n)^2$")
    axs[2].legend()
    axs[3].set_xlabel(r"$\left<z\right>/\left<\sqrt{x^2+y^2}\right>$")
    axs[3].set_ylabel(r"$\phi_0$")
    axs[3].legend()
    plt.savefig(filename[:-4]+"_slice.png")
    plt.close()
    return params[1], perr[1]


def tanalpha_q_plot():
    kcs = np.arange(0.0,1.7,0.1)
    ta = [-0.006388567029160819, 0.030965951628380804, 0.026218232784990175, 0.055499941290762016, 0.0831451415446061, 0.05713229300576396, 0.09805300889223305, 0.10824314872346479, 0.1575055186951583, 0.09059275439065695, 0.19610359491174206, 0.17714045891215047, 0.2003725358394592, 0.2254279156883871, 0.2179849393182559, 0.23040935610576754, 0.2708152045967104]

    l0=1.73
    qs = np.arctan(2*kcs/3)/(2*np.pi*(l0+1))

    fig = plt.figure(figsize=(6, 2))
    ax0 = plt.subplot2grid((1, 2), (0, 0))
    ax1 = plt.subplot2grid((1, 2), (0, 1))

    ax0.plot(kcs,ta,"o",ms=8,mfc="None")
    ax1.plot(qs,ta,"o",ms=8,mfc="None")

    ax0.set_xlabel(r"$k_c$")
    ax0.set_ylabel(r"-$\tan(\alpha)$")

    ax1.set_xlabel(r"$q = \frac{\arctan(2k_c/3)}{2\pi(l_0+1)}$")
    ax1.set_ylabel(r"-$\tan(\alpha)$")

    plt.show()
    plt.close()