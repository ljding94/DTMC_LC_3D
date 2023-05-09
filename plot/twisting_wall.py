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


def ax_pitch_phi_z_plot(filename, axcfg, axunphis, axphi0z,msize,LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4)):
    colors = ["red", "blue", "green", "tomato", "black", "purple"]

    # configuration plot
    ax_config_plot_xyz(axcfg, filename, "gray", LineWidth, pov="zx", rotxyz=(-1,0,0.5*np.pi), mesh=1, bead=0, rod=0, d=1, pwlim=np.pi / 3)

    data = np.loadtxt(filename, skiprows=6, delimiter=",", unpack=True)
    x, y, z, ux, uy, uz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    ns = np.transpose(data[17:])

    # find center (x0,y0) for the cylinderical part
    x, y, z = x - np.average(x), y - np.average(y), z - np.average(z)
    r = np.sqrt(x**2+y**2)
    zmin, zmax = np.min(z), np.max(z)
    zmin, zmax = zmin + (zmax - zmin) * z_relative_range[0], zmin + (zmax - zmin) * z_relative_range[1]
    phi = np.arctan2(y, x)
    xmax = 1.3*np.max(x)

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
            #axunphis[i].plot(phi_sort, un2_sort + i * 1, "o", markersize=msize, alpha=0.5, color=colors[i])
            axunphis[i].plot(phi_sort, un2_sort, "o", markersize=msize, alpha=0.5, color=colors[i]) # no need to lift since plotting in different axis

            # find fitting parameters
            print("phi0_bound", phi0_bound)
            para_fit = un2_fit_phi(phi_sort, un2_sort, phi0_bound)
            phi0_bound = [para_fit[2]-1.1, para_fit[2]+1.1]
            # plot fit sin curve
            #axunphis[i].plot(phi_sort, test_func_sin(phi_sort, para_fit[0], para_fit[1], para_fit[2]) + i * 1, color=colors[i], label=r"$\phi_0=%.1f$" % para_fit[2],linewidth=LineWidth)
            axunphis[i].plot(phi_sort, test_func_sin(phi_sort, para_fit[0], para_fit[1], para_fit[2]), color=colors[i], label=r"$\phi_0=%.1f$" % para_fit[2],linewidth=LineWidth) # no need to lift since plotting in different axis
            phi0.append(para_fit[2])
            print("para_fit ", para_fit)
    z_r_mean = np.array(z_mean)/np.array(R_mean)
    axphi0z.plot(z_r_mean, phi0, "o", markersize=msize, mfc="None")
    params, pcov = optimize.curve_fit(test_func_abx, z_r_mean, phi0, p0=[0, 0.1])
    perr = np.sqrt(np.diag(pcov))
    axphi0z.plot(z_r_mean, test_func_abx(z_r_mean, params[0], params[1]), "k-", label=r"$\phi_0 \propto %.2f \left<z\right>/\left<R\right>$" % params[1],linewidth=LineWidth)
    #label=r"$\phi_0 = %.1f+%.3f z$" % (params[0], params[1])


def wall_pitch_tan_data_get(date = "07Dec22"):
    kcs = np.arange(0.0,2.01,0.1)
    skips = {}
    #res by cn
    # date Dec7_2022 data
    if date == "07Dec22":
        Cn2res = [(0.07171687507132235, 0.014269194742233529), (0.1531247950928125, 0.019702645691059745), (0.21597286919392358, 0.015530465542011793), (0.12134238384877459, 0.020801319512674853), (0.2893139878856185, 0.0049807115122890485), (0.30840884105694, 0.018021602282221025), (0.3439241206076535, 0.023715579422829116), (0.3549383315516277, 0.014143177344371568), (0.3765076695151387, 0.029217039592043484), (0.48119482596213653, 0.016506420556048142), (0.6269013241538041, 0.008668224771367267), (0.6155285243118631, 0.02951455537466414), (0.6071443855656959, 0.030377613877672377), (0.6828988133368302, 0.016270617400019716), (0.7268114381338039, 0.029385986928657605), (0.807921890344102, 0.02512949370453941), (0.8116053338525985, 0.027035593644727143), (0.8498225919074379, 0.023119548395428357), (0.8019832448359602, 0.023494484374115946), (0.8416905938095114, 0.01768179574411115), (0.9213251792310383, 0.04878216325679875)]

        Cn4res = [(0.04526928659475114, 0.02450490350304673), (0.052519524383029245, 0.014833023186659607), (0.17082072926035433, 0.014317561384449597), (0.19119114042378088, 0.020541188626334566), (0.18136254424639991, 0.02977076209961673), (0.29167081623340163, 0.019907713717815103), (0.38529827242545833, 0.019665338859497696), (0.47357895083010587, 0.031999344042706206), (0.35935583051548764, 0.0035951730371498853), (0.418774417527803, 0.021991532090154482), (0.33655979350098914, 0.013728473271936882), (0.5795742763531372, 0.006700866889811913), (0.633888494763439, 0.036280314067395784), (0.7054085401292735, 0.015219188144950999), (0.7283219380194356, 0.038371045899881606), (0.7279032948946208, 0.020915320054828775), (0.8140908945784836, 0.021769667918782345), (0.9326531697794475, 0.020552343320187443), (0.9972649387446532, 0.03218122500117304), (0.9299581288732129, 0.030817098325748242), (0.9069568906852431, 0.018944770171753167)]

        Cn6res = [(-0.04535629139876542, 0.021113893611373585), (0.0859852488954683, 0.02332645619289845), (0.10046459082681683, 0.02004970392703706), (0.132936078911145, 0.029107698275464233), (0.2173186393632259, 0.00869240606734223), (0.28203634123467103, 0.022155009014534842), (0.3643730464719448, 0.004330547417869567), (0.3626924134801637, 0.018990989593462577), (0.4647405159949755, 0.009079895858261126), (0.5022692729650203, 0.033354261771154796), (0.4694281510881381, 0.014153187156549679), (0.5566644040589314, 0.027962283011405485), (0.5451418024359513, 0.033523225772062584), (0.6813221621068155, 0.014908556442604132), (0.7433417090200445, 0.04024179292816863), (0.7711906286273695, 0.014836617946756057), (0.8194082073425, 0.029523235446106667), (0.8228541736488066, 0.0198543734727229), (0.8454411583713057, 0.017040105556118172), (0.909081676268063, 0.018714415026368024), (0.907172499406838, 0.006310711006587178)]

        alpha_u_Cn2 = [1.3973005177504485, 1.4900114242662634, 1.4161562287081038, 1.3862472778494603, 1.3712475033359504, 1.2822010660616556, 1.2476022107234472, 1.2622307499337568, 1.2129765442873435, 1.1549972762064071, 1.1349823121170068, 1.133962514187405, 1.1974279096920175, 1.15748541958521, 1.117944043265488, 1.1785343323920476, 1.1453405565825097, 1.1132352415982831, 1.102632547439785, 1.082502657644467, 1.1373164644140825]
        alpha_u_Cn4 = [1.4777834484135257, 1.4551264087420979, 1.4207589949494193, 1.3716097525234503, 1.3916646625516402, 1.3072236309459868, 1.2951939281106437, 1.2112390699261928, 1.2023504703701497, 1.130248181346157, 1.0715173721614377, 1.171582291965305, 1.1596850527927227, 1.0970915891251725, 1.1536775294698973, 1.116310670673257, 1.1236886156423231, 1.1357779370983998, 1.0867051511512427, 1.1117527004756618, 1.214311917422728]
        alpha_u_Cn6 = [1.4723108811252246, 1.5041477642870125, 1.4917032885592725, 1.452903089289555, 1.4079053024350594, 1.2838809156284545, 1.2616937264686248, 1.286026793720392, 1.201167754223943, 1.1713746818378825, 1.179985018426367, 1.1796708016489579, 1.05374140829622, 1.1316281497275797, 1.1489902037650463, 1.0977612630096436, 1.1029307965274928, 1.0954793387023247, 1.1320166475595055, 1.1146776358131572, 1.120967449783173]


    elif date == "30Apr23":
        Cn2res = [(-0.1266917912297794, 0.028526144941257346), (0.03448459711010362, 0.021914611891257332), (0.14206239416067942, 0.017764021142885944), (0.22794243775152848, 0.032057432175757954), (0.1507376230559613, 0.007342720988852514), (0.3419144779124007, 0.017372991458620566), (0.32026997771143173, 0.023791505657524647), (0.38749059993344853, 0.022949144458460795), (0.5410837577031876, 0.026214282207906427), (0.5142430326951268, 0.021072181324653355), (0.49776940417296167, 0.017377034523919217), (0.6085401084765352, 0.014436498669558978), (0.7199551625493794, 0.04317562205592464), (0.6850685858717577, 0.04317093876367643), (0.6788219692688503, 0.014579482469368658), (0.7523789230604944, 0.01964163835734322), (0.7438863115931948, 0.021158990454350806), (0.7249128151697681, 0.022558553381528332), (0.08483305286422016, 0.1732957368363053), (0.9610505448774892, 0.041433492016104384), (1.0242414505969502, 0.018165144611271725)]

        alpha_u_Cn2 = [1.4536502683367474, 1.5008864602876106, 1.4298242207760699, 1.3948865379129618, 1.3403124680111602, 1.2623773236398417, 1.3052043219371414, 1.2289067903679272, 1.2359757261695994, 1.1317656543653687, 1.2055757077200768, 1.1494504378281738, 1.1760822208883062, 1.1430329316557593, 1.1063317515497106, 1.1159600769695934, 1.1151404358811252, 1.094876277560001, 1.0223460968190237, 1.1322001802555726, 1.1723219911870746]

        Cn6res = [(-0.01922937631836433, 0.01963807352899583), (0.08125109164756063, 0.013382281485955749), (0.12130483903370315, 0.034895253396119666), (0.12639319331359652, 0.008462851451805066), (0.2043066482110673, 0.0251905008878562), (0.1857073730147812, 0.02492102636366807), (0.2986785126450209, 0.011994533852154588), (0.25027504754513763, 0.011101630339116092), (0.43100642555692487, 0.017751279302732275), (0.3893497046413637, 0.015946610231412944), (0.555718005481073, 0.012377760605866935), (0.5878449085562648, 0.022549611536670504), (0.5557981903712348, 0.01961660997746607), (0.5896487272030112, 0.013473312549195553), (0.6100632889538995, 0.022016217928861423), (0.6783656456784983, 0.012951294059141121), (0.725789382184821, 0.01770455366515412), (0.7409865377932232, 0.032975443001650916), (0.7419725335284941, 0.018983540421732916), (0.9093475650376862, 0.028207926055303565), (0.9706653681279965, 0.013284244369809468)]

        alpha_u_Cn6 = [1.4714654357062324, 1.4575304430869886, 1.4067202469922706, 1.3961578240297008, 1.2924201444935788, 1.3495072401951909, 1.3330172908717124, 1.2342889372070387, 1.2467632333643517, 1.186270939204213, 1.2622193959567545, 1.1208306416293747, 1.1449681730107903, 1.1417538486145935, 1.0957334853421825, 1.1143667137000213, 1.0705347911305676, 1.0932534138488819, 1.0610275596342196, 1.138363464153026, 1.0757924546471185]

        Cn10res = [(0.5320036419942554, 0.16985677026033116), (0.7409560935086157, 0.08836030053038502), (-0.4555559374187401, 0.2161245170226461), (0.047344663804099595, 0.018522430562407353), (0.1252471388348498, 0.013226304820043448), (0.14305276311087384, 0.026718288029681694), (0.21325604437544712, 0.01325945077681793), (0.2620877440107482, 0.013314947076784133), (0.33569171691069577, 0.020484705398043302), (0.4575902439057407, 0.012276430931733936), (0.3907546063390977, 0.04296482633730491), (0.5922038895427746, 0.016169677009297785), (0.5387520536725244, 0.022684990772079118), (0.5659818661793371, 0.026805302422244544), (0.6281478094197216, 0.025035669793035724), (0.7039746814278789, 0.012660329725113953), (0.8272223033389206, 0.017560732609299273), (0.7879447563994834, 0.012784495007270401), (0.8099440155707626, 0.020722140839961364), (0.7953545915942407, 0.012351842158211846), (0.8300654466064288, 0.038448685070608196)]

        alpha_u_Cn10 = [-100, -100, -100, 1.412646452801134, 1.3641394550966177, 1.343791066633191, 1.2346498382504412, 1.2757218181447256, 1.1515589840806653, 1.1396427962946054, 1.139608973026255, 1.1521907009625163, 1.1748493961842763, 1.0973241955431385, 1.1679076273294344, 1.1115510242401012, 1.1511004505454516, 1.080418521356136, 1.0755006110710026, 1.111711364873674, 1.0996561668047167]

        skips = {"Cn2":2,"Cn6":1,"Cn10":4}


    return [kcs, Cn2res,Cn6res,Cn10res,alpha_u_Cn2,alpha_u_Cn6,alpha_u_Cn10,skips]



def wall_pitch_q_plot(LineWidth, FontSize, LabelSize):
    print("👌")
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.3))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((15, 2), (0, 0), rowspan=6)
    axunphi0 = plt.subplot2grid((15, 2), (6, 0), rowspan=9)
    divider = make_axes_locatable(axunphi0)
    axunphi1 = divider.append_axes("top",size="100%",pad=0.03)
    axunphi2 = divider.append_axes("top",size="100%",pad=0.03)
    axunphi3 = divider.append_axes("top",size="100%",pad=0.03)
    axunphi4 = divider.append_axes("top",size="100%",pad=0.03)
    axunphis = [axunphi0,axunphi1,axunphi2,axunphi3,axunphi4]
    #for i in range(5):
        #axunphis.append()

    axphi0z = plt.subplot2grid((15, 2), (0, 1), rowspan=5)
    axtankc = plt.subplot2grid((15, 2), (5, 1), rowspan=5)
    axthekc = plt.subplot2grid((15, 2), (10, 1), rowspan=5, sharex = axtankc)

    msize = 4

    # single tube plot
    #fname = "../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn4.0.csv"
    #fname = "../data/Ne2/Dec7_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn4.0_id0.csv"
    #fname = "../data/Ne2/Dec7_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn6.0_id0.csv"
    fname = "../data/Ne2/Apr30_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.6_Cn10.0_id0.csv"

    ax_pitch_phi_z_plot(fname, axcfg, axunphis, axphi0z, msize, LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4))
    axcfg.margins(y=0)
    # config plot
    x1, y1 = 0.7, 1.2+0.075/6*5
    axcfg.text(x1, y1, r"(a)", fontsize=FontSize, transform=axunphis[4].transAxes)

    # print(filename, params[1])
    # axs[2].scatter(phi,x)
    for i in range(5):
        #axunphis[i].yaxis.set_label_position("left")
        axunphis[i].tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=False, labelleft=True, labelsize=LabelSize)
        #axunphis[i].set_xlim(-2.7, 2.7)
        axunphis[i].yaxis.set_major_locator(MultipleLocator(0.5))
        axunphis[i].yaxis.set_minor_locator(MultipleLocator(0.25))
        axunphis[i].set_ylim(-0.1, 1.1)
    axunphi0.set_xlabel(r"$\phi = \arctan(y/x)$",fontsize=FontSize)
    axunphi0.xaxis.set_major_locator(MultipleLocator(2))
    axunphi0.xaxis.set_minor_locator(MultipleLocator(1))
    axunphi0.tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=True, labelleft=True, labelsize=LabelSize)
    axunphis[2].set_ylabel(r"$(\vu{u}\cdot \vu{n})^2$",fontsize=FontSize)

    x1, y1 = 0.7, 0.057*5/9*5
    axunphi0.text(x1, y1, r"(b)", fontsize=FontSize, transform=axunphis[0].transAxes)


    axphi0z.set_xlabel(r"$\left<z\right>/\left<R\right>$",fontsize=FontSize)
    axphi0z.set_ylabel(r"$\phi_0$",fontsize=FontSize)
    axphi0z.yaxis.set_label_position("left")
    #axphi0z.yaxis.tick_right()
    axphi0z.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=True, labelsize=LabelSize)

    axphi0z.xaxis.set_major_locator(MultipleLocator(1))
    axphi0z.xaxis.set_minor_locator(MultipleLocator(0.5))
    axphi0z.set_xlim(-2.15, 2.15)
    axphi0z.set_ylim(0.7, 1.65)
    axphi0z.yaxis.set_major_locator(MultipleLocator(0.2))
    axphi0z.yaxis.set_minor_locator(MultipleLocator(0.1))

    axphi0z.legend(loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)

    x1, y1 = 0.75, 0.075
    axphi0z.text(x1, y1, r"(c)", fontsize=FontSize, transform=axphi0z.transAxes)

    # no longer needed, since R already considered in fitting slope
    kcs, Cn2res,Cn6res,Cn10res,alpha_u_Cn2,alpha_u_Cn6,alpha_u_Cn10, skips = wall_pitch_tan_data_get(date = "30Apr23")
    R = 1
    Cn2res = R*np.array(Cn2res)
    Cn2ta, Cn2taerr = Cn2res[:,0], Cn2res[:,1]
    Cn6res = R*np.array(Cn6res)
    Cn6ta, Cn6taerr = Cn6res[:,0], Cn6res[:,1]
    Cn10res = R*np.array(Cn10res)
    Cn10ta, Cn10taerr = Cn10res[:,0], Cn10res[:,1]

    #kcs = np.delete(kcs,9)
    #ta = np.delete(ta,9)
    #taerr = np.delete(taerr,9)

    #print(ta)
    #print(taerr)

    #axtankc.errorbar(kcs,ta,yerr = taerr, color = "k", marker = "o", linestyle="None" ,ms=msize,mfc="None")
    n0Cn2,n0Cn6,n0Cn10 = skips["Cn2"],skips["Cn6"],skips["Cn10"]
    axtankc.errorbar(np.delete(kcs[n0Cn2:],-3),np.delete(Cn2ta[n0Cn2:],-3),yerr = np.delete(Cn2taerr[n0Cn2:],-3), color = "b", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=2")
    axtankc.errorbar(kcs[n0Cn6:],Cn6ta[n0Cn6:],yerr = Cn6taerr[n0Cn6:], color = "r", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=6")
    axtankc.errorbar(kcs[n0Cn10:],Cn10ta[n0Cn10:],yerr = Cn10taerr[n0Cn10:], color = "k", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=10")
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

    axthekc.plot(np.delete(kcs[n0Cn2:],-3),np.delete(alpha_u_Cn2[n0Cn2:],-3) - np.arctan(np.delete(Cn2ta[n0Cn2:],-3)), color = "b", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=2")
    axthekc.plot(kcs[n0Cn6:],alpha_u_Cn6[n0Cn6:] - np.arctan(Cn6ta[n0Cn6:]), color = "r", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=6")
    axthekc.plot(kcs[n0Cn10:],alpha_u_Cn10[n0Cn10:] - np.arctan(Cn10ta[n0Cn10:]), color = "k", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=10")
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

    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.4,hspace=0.3,left=0.12,right=0.97,bottom=0.08,top=0.98)
    plt.savefig("figures/wall_pitch.pdf", format="pdf")
    plt.close()
