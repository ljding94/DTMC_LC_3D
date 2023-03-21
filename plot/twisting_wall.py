import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from scipy import optimize


def twist_q_data_get():
    # foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    # foldername = "../data/Ne2/Mar_2022/Mar23_2022"
    foldername = "../data/Ne2/May12_2022"
    lfs = [15, 25, 35]
    Kd = 4.0
    Cn = 4.0
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue", "purple"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(lfs)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_qs_Cn%.1f_ana.csv" % (lfs[i], Kd, Cn)
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    qs, uc_aves, uc_errs, uz2_aves, uz2_errs = datas[0], datas[28], datas[30], datas[34], datas[36]
    labels = list(map(str, lfs))
    legendtitle = r"$l_f$"
    return [qs, uc_aves, uc_errs, uz2_aves, uz2_errs, labels, colors, markers, legendtitle]


def twist_q_config_data_get():
    # foldername = "../data/Ne2/Mar_2022/Mar23_2022"
    foldername = "../data/Ne2/May12_2022"
    lf = 25.0
    Kd = 4.0
    Cn = 4.0
    lfs = [25, 25, 25, 35, 35, 35]
    qs = [0, 1, 2.5, 0, 1, 2.5]
    fnames, povs, rotxyzs, xysfts = [], [], [], []
    for i in range(len(qs)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f.csv" % (lfs[i], Kd, qs[i], Cn))
        povs.append("zx")
        # rotxyzs.append([np.pi/3,0,np.pi/2])
        rotxyzs.append([np.pi / 3, 0, 0])
    rotxyzs[3] = [np.pi / 2 - 0.1, 0, 0]
    xysfts = [[0, 0], [0, 13], [0, 26], [35, 0], [35, 13], [35, 26]]
    return [lfs, qs, fnames, povs, rotxyzs, xysfts]


def twist_q_plot(LineWidth, FontSize, LabelSize):
    print("👌")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 2, 246 / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((2, 10), (0, 0), colspan=7, rowspan=2)
    # axcfgt = plt.subplot2grid((2, 5), (1, 3),colspan=2)
    axuc = plt.subplot2grid((2, 10), (0, 7), colspan=3)
    axuz = plt.subplot2grid((2, 10), (1, 7), colspan=3, sharex=axuc)
    msize = 4

    """
    # twist bond color
    for i in range(len(qs)):
        pass
        ax_config_plot_xyz(axcfgt, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=0,yshift=15*i, mesh=1, bead=0,rod=0,twistcolor=1)
        axcfgt.text(-15, 15*i-7,r"$k_c=%.0f$"%qs[i],fontsize=FontSize)
    cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.25*np.pi), cmap=cm.get_cmap("ocean")),ax=axcfgt,ticks=[0,0.25*np.pi,0.5*np.pi],orientation="horizontal")
    cbar.ax.set_xticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"],fontsize=FontSize)
    cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    cbar.ax.set_title(r"$T_s$",fontsize=FontSize)


    axcfgt.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    x1, y1 = 0.8, -0.0
    axcfgt.text(x1,y1, r"(d)", fontsize=FontSize,transform=axcfgt.transAxes)
    """

    # configuration plot
    lfs, qs, fnames, povs, rotxyzs, xysfts = twist_q_config_data_get()
    for i in range(len(qs)):
        pass
        print(xysfts[i])
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i], xshift=xysfts[i][0], yshift=xysfts[i][1], mesh=1, bead=0, rod=0, d=1, pwlim=np.pi / 3)
        axcfg.text(xysfts[i][0] - lfs[i] / 2, xysfts[i][1] - 6, r"$k_c=%.1f,l_f=%.0f$" % (qs[i], lfs[i]), fontsize=FontSize)

        # ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, xshift=10,yshift=25*i-5,zslice=(8,12), mesh=1, bead=0,rod=0,pwlim=0.8, d=1)
    #axcfg.text(xysfts[2][0] - 5, xysfts[2][1] + lfs[i] / 2 - 1, r"$k_c=%.1f$" % qs[-1], fontsize=FontSize)
    axcfg.tick_params(which="both", direction="in", bottom="off", top="off", right="off", left="off", labelbottom=False, labelleft=False, labelsize=LabelSize)
    x1, y1 = -0.5, 0.1
    axcfg.text(x1, y1, r"(a)", fontsize=FontSize, transform=axuz.transAxes)
    # axcfg.text(x1,y1, r"(c)", fontsize=FontSize,transform=axcfg.transAxes)


    # o vs q plot

    qs, uc_aves, uc_errs, uz2_aves, uz2_errs, labels, colors, markers, legendtitle = twist_q_data_get()
    nf = 30
    n = 2
    for i in range(len(qs)):
        axuc.errorbar(qs[i][:nf:n], uc_aves[i][:nf:n], uc_errs[i][:nf:n], ls="None", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])

    axuc.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=LabelSize)
    # axuc.set_ylabel(r"$\left<(\vu{u}_i\cross\vu{u}_j)\cdot \vu{r}_{ij} (\vu{u}_i\cdot\vu{u}_j)\right>_{(i,j)}$", fontsize=FontSize)
    axuc.set_ylabel(r"$T_s$", fontsize=FontSize)
    # axuc.set_ylim(0.0,0.22)
    axuc.xaxis.set_major_locator(MultipleLocator(1.0))
    axuc.xaxis.set_minor_locator(MultipleLocator(0.5))
    axuc.yaxis.set_major_locator(MultipleLocator(0.05))
    axuc.yaxis.set_minor_locator(MultipleLocator(0.025))
    # axuc.set_xlabel(r"$k_c$",fontsize=FontSize)
    axuc.legend(loc="upper left", title=legendtitle, ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    x1, y1 = 0.8, 0.1
    axuc.text(x1, y1, r"(b)", fontsize=FontSize, transform=axuc.transAxes)

    for i in range(len(qs)):
        axuz.errorbar(qs[i][:nf:n], uz2_aves[i][:nf:n], uz2_errs[i][:nf:n], ls="None", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])

    axuz.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=LabelSize)
    axuz.set_ylabel(r"$(\vu{u}\cdot \vu{z})^2$", fontsize=FontSize)
    # axuz.set_ylim(0.0,0.2)
    axuz.xaxis.set_major_locator(MultipleLocator(0.5))
    axuz.xaxis.set_minor_locator(MultipleLocator(0.25))
    axuz.set_xlim(-0.2, 3.2)
    axuz.yaxis.set_major_locator(MultipleLocator(0.05))
    axuz.yaxis.set_minor_locator(MultipleLocator(0.025))
    axuz.set_xlabel(r"$k_c$", fontsize=FontSize)
    axuz.legend(loc="upper left", title=legendtitle, ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    x1, y1 = 0.8, 0.1
    axuz.text(x1, y1, r"(c)", fontsize=FontSize, transform=axuz.transAxes)


    # fig.tight_layout(pad=0.1)
    plt.tight_layout(pad=0.1)
    plt.savefig("figures/twisting_wall.pdf", format="pdf")

    # extra plot just to see the trend agains 1/2arctan(kc/3)
    """
    plt.figure()
    for i in range(len(qs)):
        axuz.errorbar(np.arctan(qs[i][:nf:n]/3)/2,np.sqrt(uz2_aves[i][:nf:n]),uz2_errs[i][:nf:n], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    plt.show()
    plt.close()
    """


def test_func_sin(phi, a, b, c):
    return a + 0.5 * np.sin(2 * (phi - c))


def test_func_abx(x, a, b):
    return a + b * np.array(x)


def un2_fit_phi(phi, un2, phi0_bound):
    params, params_covariance = optimize.curve_fit(test_func_sin, phi, un2, p0=[0.5, 0.5, np.average(phi0_bound)], bounds=((0.45, 0.45, phi0_bound[0]), (0.55, 0.55, phi0_bound[1])))
    print("params: ", params)
    print("params_covariance: ", params_covariance)
    return params


def ax_pitch_phi_z_plot(filename, axcfg, axunphi, axphi0z,msize,LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4)):
    colors = ["red", "blue", "green", "tomato", "black", "purple"]

    # configuration plot
    ax_config_plot_xyz(axcfg, filename, "gray", LineWidth, pov="zx", rotxyz=(0,0,0.5*np.pi), mesh=1, bead=0, rod=0, d=1, pwlim=np.pi / 3)

    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x, y, z, ux, uy, uz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    ns = np.transpose(data[17:])

    # find center (x0,y0) for the cylinderical part
    x, y, z = x - np.average(x), y - np.average(y), z - np.average(z)
    r = np.sqrt(x**2+y**2)
    zmin, zmax = np.min(z), np.max(z)
    zmin, zmax = zmin + (zmax - zmin) * z_relative_range[0], zmin + (zmax - zmin) * z_relative_range[1]
    phi = np.arctan2(y, x)
    xmax = 1.1*np.max(x)

    z_mean = []
    R_mean = []
    phi0 = []
    for i in range(nbin):
        zl, zr = zmin + (zmax - zmin) * i / nbin, zmin + (zmax - zmin) * (i + 1) / nbin
        select = np.logical_and(z > zl, z <= zr)
        z_mean.append(np.average(z[select]))
        R_mean.append(np.average(r[select]))

        # add range indicator to config plot
        axcfg.plot([xmax,xmax],[zl,zr],color=colors[i],linewidth=3,solid_capstyle="butt")
        axcfg.plot([-xmax,xmax],[zl,zl],"k--",alpha=0.2,linewidth=1)
        axcfg.plot([-xmax,xmax],[zr,zr],"k--",alpha=0.2,linewidth=1)

        if i == 0:
            phi0_bound = [-np.pi/2, np.pi/2]
        if 1:
            un2_sort = un2[select][phi[select].argsort()]
            phi_sort = np.sort(phi[select])
            # plot points
            axunphi.plot(phi_sort, un2_sort + i * 1, "o", markersize=msize, alpha=0.5, color=colors[i])

            # find fitting parameters
            print("phi0_bound", phi0_bound)
            para_fit = un2_fit_phi(phi_sort, un2_sort, phi0_bound)
            phi0_bound = [para_fit[2]-1.1, para_fit[2]+1.1]
            # plot fit sin curve
            axunphi.plot(phi_sort, test_func_sin(phi_sort, para_fit[0], para_fit[1], para_fit[2]) + i * 1, color=colors[i], label=r"$\phi_0=%.1f$" % para_fit[2],linewidth=LineWidth)
            phi0.append(para_fit[2])
            print("para_fit ", para_fit)
    z_r_mean = np.array(z_mean)/np.array(R_mean)
    axphi0z.plot(z_r_mean, phi0, "o", markersize=msize, mfc="None")
    params, pcov = optimize.curve_fit(test_func_abx, z_r_mean, phi0, p0=[0, 0.1])
    perr = np.sqrt(np.diag(pcov))
    axphi0z.plot(z_r_mean, test_func_abx(z_r_mean, params[0], params[1]), "k-", label=r"$\phi_0 \propto %.2f \left<z\right>/\left<R\right>$" % params[1],linewidth=LineWidth)
    #label=r"$\phi_0 = %.1f+%.3f z$" % (params[0], params[1])



def wall_pitch_tan_data_get():
    kcs = np.arange(0.0,2.01,0.1)
    #res by cn
    Cn2res = [(0.07171687507132235, 0.014269194742233529), (0.1531247950928125, 0.019702645691059745), (0.21597286919392358, 0.015530465542011793), (0.12134238384877459, 0.020801319512674853), (0.2893139878856185, 0.0049807115122890485), (0.30840884105694, 0.018021602282221025), (0.3439241206076535, 0.023715579422829116), (0.3549383315516277, 0.014143177344371568), (0.3765076695151387, 0.029217039592043484), (0.48119482596213653, 0.016506420556048142), (0.6269013241538041, 0.008668224771367267), (0.6155285243118631, 0.02951455537466414), (0.6071443855656959, 0.030377613877672377), (0.6828988133368302, 0.016270617400019716), (0.7268114381338039, 0.029385986928657605), (0.807921890344102, 0.02512949370453941), (0.8116053338525985, 0.027035593644727143), (0.8498225919074379, 0.023119548395428357), (0.8019832448359602, 0.023494484374115946), (0.8416905938095114, 0.01768179574411115), (0.9213251792310383, 0.04878216325679875)]

    Cn4res = [(0.04526928659475114, 0.02450490350304673), (0.052519524383029245, 0.014833023186659607), (0.17082072926035433, 0.014317561384449597), (0.19119114042378088, 0.020541188626334566), (0.18136254424639991, 0.02977076209961673), (0.29167081623340163, 0.019907713717815103), (0.38529827242545833, 0.019665338859497696), (0.47357895083010587, 0.031999344042706206), (0.35935583051548764, 0.0035951730371498853), (0.418774417527803, 0.021991532090154482), (0.33655979350098914, 0.013728473271936882), (0.5795742763531372, 0.006700866889811913), (0.633888494763439, 0.036280314067395784), (0.7054085401292735, 0.015219188144950999), (0.7283219380194356, 0.038371045899881606), (0.7279032948946208, 0.020915320054828775), (0.8140908945784836, 0.021769667918782345), (0.9326531697794475, 0.020552343320187443), (0.9972649387446532, 0.03218122500117304), (0.9299581288732129, 0.030817098325748242), (0.9069568906852431, 0.018944770171753167)]

    Cn6res = [(-0.04535629139876542, 0.021113893611373585), (0.0859852488954683, 0.02332645619289845), (0.10046459082681683, 0.02004970392703706), (0.132936078911145, 0.029107698275464233), (0.2173186393632259, 0.00869240606734223), (0.28203634123467103, 0.022155009014534842), (0.3643730464719448, 0.004330547417869567), (0.3626924134801637, 0.018990989593462577), (0.4647405159949755, 0.009079895858261126), (0.5022692729650203, 0.033354261771154796), (0.4694281510881381, 0.014153187156549679), (0.5566644040589314, 0.027962283011405485), (0.5451418024359513, 0.033523225772062584), (0.6813221621068155, 0.014908556442604132), (0.7433417090200445, 0.04024179292816863), (0.7711906286273695, 0.014836617946756057), (0.8194082073425, 0.029523235446106667), (0.8228541736488066, 0.0198543734727229), (0.8454411583713057, 0.017040105556118172), (0.909081676268063, 0.018714415026368024), (0.907172499406838, 0.006310711006587178)]

    alpha_u_Cn2 = [1.3973005177504485, 1.4900114242662634, 1.4161562287081038, 1.3862472778494603, 1.3712475033359504, 1.2822010660616556, 1.2476022107234472, 1.2622307499337568, 1.2129765442873435, 1.1549972762064071, 1.1349823121170068, 1.133962514187405, 1.1974279096920175, 1.15748541958521, 1.117944043265488, 1.1785343323920476, 1.1453405565825097, 1.1132352415982831, 1.102632547439785, 1.082502657644467, 1.1373164644140825]
    alpha_u_Cn4 = [1.4777834484135257, 1.4551264087420979, 1.4207589949494193, 1.3716097525234503, 1.3916646625516402, 1.3072236309459868, 1.2951939281106437, 1.2112390699261928, 1.2023504703701497, 1.130248181346157, 1.0715173721614377, 1.171582291965305, 1.1596850527927227, 1.0970915891251725, 1.1536775294698973, 1.116310670673257, 1.1236886156423231, 1.1357779370983998, 1.0867051511512427, 1.1117527004756618, 1.214311917422728]
    alpha_u_Cn6 = [1.4723108811252246, 1.5041477642870125, 1.4917032885592725, 1.452903089289555, 1.4079053024350594, 1.2838809156284545, 1.2616937264686248, 1.286026793720392, 1.201167754223943, 1.1713746818378825, 1.179985018426367, 1.1796708016489579, 1.05374140829622, 1.1316281497275797, 1.1489902037650463, 1.0977612630096436, 1.1029307965274928, 1.0954793387023247, 1.1320166475595055, 1.1146776358131572, 1.120967449783173]

    return [kcs, Cn2res,Cn4res,Cn6res,alpha_u_Cn2,alpha_u_Cn4,alpha_u_Cn6]



def wall_pitch_q_plot(LineWidth, FontSize, LabelSize):
    print("👌")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((6, 2), (0, 0), rowspan=3)
    axunphi = plt.subplot2grid((6, 2), (3, 0), rowspan=3)
    axphi0z = plt.subplot2grid((6, 2), (0, 1), rowspan=2)
    axtankc = plt.subplot2grid((6, 2), (2, 1), rowspan=2)
    axthekc = plt.subplot2grid((6, 2), (4, 1), rowspan=2, sharex = axtankc)
    msize = 4

    # single tube plot
    #fname = "../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn4.0.csv"
    fname = "../data/Ne2/Dec7_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn4.0_id0.csv"

    ax_pitch_phi_z_plot(fname, axcfg, axunphi, axphi0z, msize, LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4))
    axcfg.margins(y=0)
    # config plot
    x1, y1 = 0.75, 1.075
    axcfg.text(x1, y1, r"(a)", fontsize=FontSize, transform=axunphi.transAxes)

    # print(filename, params[1])
    # axs[2].scatter(phi,x)
    axunphi.set_xlabel(r"$\phi = \arctan(y/x)$",fontsize=FontSize)
    axunphi.set_ylabel(r"$(\vu{u}\cdot \vu{n})^2$",fontsize=FontSize)
    axunphi.yaxis.set_label_position("left")
    #axunphi.yaxis.tick_right()
    axunphi.tick_params(which="both", direction="in", bottom="on", top="off", right="off", left="oon", labelbottom=True, labelleft=True, labelsize=LabelSize)
    axunphi.xaxis.set_major_locator(MultipleLocator(2))
    axunphi.xaxis.set_minor_locator(MultipleLocator(1))
    #axunphi.set_xlim(-2.7, 2.7)
    axunphi.yaxis.set_major_locator(MultipleLocator(1))
    axunphi.yaxis.set_minor_locator(MultipleLocator(0.5))
    axunphi.set_ylim(-0.4, 5.2)

    x1, y1 = 0.75, 0.075
    axunphi.text(x1, y1, r"(b)", fontsize=FontSize, transform=axunphi.transAxes)


    axphi0z.set_xlabel(r"$\left<z\right>/\left<R\right>$",fontsize=FontSize)
    axphi0z.set_ylabel(r"$\phi_0$",fontsize=FontSize)
    axphi0z.yaxis.set_label_position("left")
    #axphi0z.yaxis.tick_right()
    axphi0z.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=True, labelsize=LabelSize)

    axphi0z.xaxis.set_major_locator(MultipleLocator(1))
    axphi0z.xaxis.set_minor_locator(MultipleLocator(0.5))
    axphi0z.set_xlim(-2.15, 2.15)
    axphi0z.set_ylim(-1.50, 1.6)
    axphi0z.yaxis.set_major_locator(MultipleLocator(0.5))
    axphi0z.yaxis.set_minor_locator(MultipleLocator(0.25))

    axphi0z.legend(loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)

    x1, y1 = 0.75, 0.075
    axphi0z.text(x1, y1, r"(c)", fontsize=FontSize, transform=axphi0z.transAxes)

    # no longer needed, since R already considered in fitting slope
    kcs, Cn2res,Cn4res,Cn6res,alpha_u_Cn2,alpha_u_Cn4,alpha_u_Cn6 = wall_pitch_tan_data_get()
    R = 1
    Cn2res = R*np.array(Cn2res)
    Cn2ta, Cn2taerr = Cn2res[:,0], Cn2res[:,1]
    Cn4res = R*np.array(Cn4res)
    Cn4ta, Cn4taerr = Cn4res[:,0], Cn4res[:,1]
    Cn6res = R*np.array(Cn6res)
    Cn6ta, Cn6taerr = Cn6res[:,0], Cn6res[:,1]

    #kcs = np.delete(kcs,9)
    #ta = np.delete(ta,9)
    #taerr = np.delete(taerr,9)

    #print(ta)
    #print(taerr)

    #axtankc.errorbar(kcs,ta,yerr = taerr, color = "k", marker = "o", linestyle="None" ,ms=msize,mfc="None")
    axtankc.errorbar(kcs,Cn2ta,yerr = Cn2taerr, color = "b", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=2")
    axtankc.errorbar(kcs,Cn4ta,yerr = Cn4taerr, color = "r", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=4")
    axtankc.errorbar(kcs,Cn6ta,yerr = Cn6taerr, color = "k", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=6")
    axtankc.set_xlabel(r"$k_c$",fontsize=FontSize)
    axtankc.set_ylabel(r"$\tan\alpha$",fontsize=FontSize)
    axtankc.yaxis.set_label_position("left")
    #axtankc.yaxis.tick_right()
    axtankc.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="on", labelbottom=True, labelleft=True, labelsize=LabelSize)

    axtankc.xaxis.set_major_locator(MultipleLocator(0.5))
    axtankc.xaxis.set_minor_locator(MultipleLocator(0.25))
    #axtankc.set_xlim(-2.7, 2.7)
    axtankc.yaxis.set_major_locator(MultipleLocator(0.5))
    axtankc.yaxis.set_minor_locator(MultipleLocator(0.1))
    axtankc.set_ylim(-0.01,1.1)
    axtankc.legend(loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)

    x1, y1 = 0.75, 0.075
    axtankc.text(x1, y1, r"(d)", fontsize=FontSize, transform=axtankc.transAxes)

    axthekc.plot(kcs,alpha_u_Cn2 - np.arctan(Cn2ta), color = "b", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=2")
    axthekc.plot(kcs,alpha_u_Cn4 - np.arctan(Cn4ta), color = "r", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=4")
    axthekc.plot(kcs,alpha_u_Cn4 - np.arctan(Cn4ta), color = "k", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=6")
    axthekc.set_xlabel(r"$k_c$",fontsize=FontSize)
    axthekc.set_ylabel(r"$\theta = \alpha_u - \alpha$",fontsize=FontSize)
    axthekc.yaxis.set_label_position("left")
    #axthekc.yaxis.tick_right()
    axthekc.tick_params(which="both", direction="in", bottom="on", top="off", right="off", left="on", labelbottom=True, labelleft=True, labelsize=LabelSize)

    axthekc.xaxis.set_major_locator(MultipleLocator(0.5))
    axthekc.xaxis.set_minor_locator(MultipleLocator(0.25))
    axthekc.yaxis.set_major_locator(MultipleLocator(0.5))
    axthekc.yaxis.set_minor_locator(MultipleLocator(0.1))
    axthekc.set_ylim(-0.01,1.51)
    axthekc.legend(loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    axthekc.text(0.75, 0.075, r"(e)", fontsize=FontSize, transform=axthekc.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/wall_pitch.pdf", format="pdf")
    plt.close()
