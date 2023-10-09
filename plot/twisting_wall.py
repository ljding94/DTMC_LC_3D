import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from scipy import optimize


def twist_q_data_get():
    # foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    # foldername = "../data/Ne2/Mar_2022/Mar23_2022"
    # foldername = "../data/Ne2/data_2022/May12_2022"
    foldername = "../data/Ne2/Jun29_2023"
    lfs = [15, 25, 35]
    Kd = 4.0
    Cn = 4.0
    datas, labels, colors, markers = [], [], [], []
    #colors = ["red", "green", "blue", "royalblue", "purple"]
    colors = ["blue", "orange", "purple", "red"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(lfs)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_qs_Cn%.1f_id0_ana.csv" % (lfs[i], Kd, Cn)
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    qs, uc_aves, uc_errs, uz2_aves, uz2_errs = datas[0], datas[28], datas[30], datas[34], datas[36]
    labels = list(map(str, lfs))
    legendtitle = r"$l_f$"
    return [qs, uc_aves, uc_errs, uz2_aves, uz2_errs, labels, colors, markers, legendtitle]


def twist_q_config_data_get():
    # foldername = "../data/Ne2/Mar_2022/Mar23_2022"
    #foldername = "../data/Ne2/data_2022/May12_2022"
    foldername = "../data/Ne2/Jun29_2023"
    lf = 25.0
    Kd = 4.0
    Cn = 4.0
    lfs = [25, 25, 25, 35, 35, 35]
    qs = [0, 1, 2.5, 0, 1, 2.5]
    fnames, povs, rotxyzs, xysfts = [], [], [], []
    for i in range(len(qs)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f_id0.csv" % (lfs[i], Kd, qs[i], Cn))
        povs.append("zx")
        # rotxyzs.append([np.pi/3,0,np.pi/2])
        rotxyzs.append([np.pi / 3, 0, 0])
    rotxyzs[3] = [np.pi / 2 - 0.1, 0, 0]
    xysfts = [[0, 0], [0, 13], [0, 26], [35, 0], [35, 13], [35, 26]]
    return [lfs, qs, fnames, povs, rotxyzs, xysfts]


def twist_q_plot(LineWidth, FontSize, LabelSize):
    print("👌 twisting wall 2 to 3 walls versus lfs")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 2, 246 / ppi * 0.9))
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
    axuz.set_ylabel(r"$\left<(\vu{u}\cdot \vu{z})^2\right>$", fontsize=FontSize)
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

def test_func_exp_sin(phi,lamphi,phi0):
    return (np.exp(1/lamphi)-np.exp(np.sin(2*(phi-phi0))/lamphi))/(np.exp(1/lamphi)-np.exp(-1/lamphi))

def test_func_exp_sin_m3(phi,lamphi,phi0):
    return (np.exp(1/lamphi)-np.exp(np.sin(3*(phi-phi0))/lamphi))/(np.exp(1/lamphi)-np.exp(-1/lamphi))

def un2_fit_phi(phi, un2, phi0_bound,m=2):
    #params, params_covariance = optimize.curve_fit(test_func_sin, phi, un2, p0=[0.5, 0.5, np.average(phi0_bound)], bounds=((0.45, 0.45, phi0_bound[0]), (0.55, 0.55, phi0_bound[1])))
    if m==3:
        params, params_covariance = optimize.curve_fit(test_func_exp_sin_m3, phi, un2, p0=[5, np.average(phi0_bound)],bounds=((0.02,phi0_bound[0]),(20,phi0_bound[1])))
    else:
        params, params_covariance = optimize.curve_fit(test_func_exp_sin, phi, un2, p0=[5, np.average(phi0_bound)],bounds=((0.02,phi0_bound[0]),(20,phi0_bound[1])))
    print("params: ", params)
    print("params_covariance: ", params_covariance)
    return params


def ax_pitch_phi_z_plot(filename, axcfg, axunphis, axphi0z,msize,LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4),m=2):
    colors = ["red", "blue", "purple", "tomato", "black"]

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
            phi0_bound = [-3*np.pi/4,np.pi/4]
        if 1:
            un2_sort = un2[select][phi[select].argsort()]
            phi_sort = np.sort(phi[select])
            # plot points
            #axunphis[i].plot(phi_sort, un2_sort + i * 1, "o", markersize=msize, alpha=0.5, color=colors[i])
            axunphis[i].plot(phi_sort, un2_sort, "o", markersize=msize, alpha=0.5, color=colors[i]) # no need to lift since plotting in different axis

            # find fitting parameters
            print("phi0_bound", phi0_bound)
            para_fit = un2_fit_phi(phi_sort, un2_sort, phi0_bound,m)
            phi0_bound = [para_fit[1]-1, para_fit[1]+1]

            # plot fit exp sin curve
            phi_p = np.linspace(-np.pi,np.pi,100)
            #axunphis[i].plot(phi_sort, test_func_sin(phi_sort, para_fit[0], para_fit[1], para_fit[2]) + i * 1, color=colors[i], label=r"$\phi_0=%.1f$" % para_fit[2],linewidth=LineWidth)
            if m==3:
                axunphis[i].plot(phi_p, test_func_exp_sin_m3(phi_p, para_fit[0], para_fit[1]), color=colors[i], label=r"$\phi_0=%.1f$" % para_fit[1],linewidth=LineWidth) # no need to lift since plotting in different axis
            else:
                axunphis[i].plot(phi_p, test_func_exp_sin(phi_p, para_fit[0], para_fit[1]), color=colors[i], label=r"$\phi_0=%.1f$" % para_fit[1],linewidth=LineWidth) # no need to lift since plotting in different axis
            phi0.append(para_fit[1])
            print("para_fit ", para_fit)

    if(axphi0z):
        z_r_mean = np.array(z_mean)/np.array(R_mean)
        axphi0z.plot(z_r_mean, phi0, "o", markersize=msize, mfc="None")
        #params, pcov = optimize.curve_fit(test_func_abx, z_r_mean, phi0, p0=[0, 0.1])
        params, pcov = optimize.curve_fit(test_func_abx, z_r_mean, phi0, p0=[0, 1])
        perr = np.sqrt(np.diag(pcov))
        axphi0z.plot(z_r_mean, test_func_abx(z_r_mean, params[0], params[1]), "k-", label=r"$\phi_0 \propto %.2f \left<z\right>/\left<R\right>$" % params[1],linewidth=LineWidth)
        #label=r"$\phi_0 = %.1f+%.3f z$" % (params[0], params[1])


def wall_pitch_tan_data_get(date = "07Dec22"):

    skips = {}
    #res by cn
    # date Dec7_2022 data
    if date == "07Dec22":
        kcs = np.arange(0.0,2.01,0.1)
        Cn2res = [(0.07171687507132235, 0.014269194742233529), (0.1531247950928125, 0.019702645691059745), (0.21597286919392358, 0.015530465542011793), (0.12134238384877459, 0.020801319512674853), (0.2893139878856185, 0.0049807115122890485), (0.30840884105694, 0.018021602282221025), (0.3439241206076535, 0.023715579422829116), (0.3549383315516277, 0.014143177344371568), (0.3765076695151387, 0.029217039592043484), (0.48119482596213653, 0.016506420556048142), (0.6269013241538041, 0.008668224771367267), (0.6155285243118631, 0.02951455537466414), (0.6071443855656959, 0.030377613877672377), (0.6828988133368302, 0.016270617400019716), (0.7268114381338039, 0.029385986928657605), (0.807921890344102, 0.02512949370453941), (0.8116053338525985, 0.027035593644727143), (0.8498225919074379, 0.023119548395428357), (0.8019832448359602, 0.023494484374115946), (0.8416905938095114, 0.01768179574411115), (0.9213251792310383, 0.04878216325679875)]

        Cn4res = [(0.04526928659475114, 0.02450490350304673), (0.052519524383029245, 0.014833023186659607), (0.17082072926035433, 0.014317561384449597), (0.19119114042378088, 0.020541188626334566), (0.18136254424639991, 0.02977076209961673), (0.29167081623340163, 0.019907713717815103), (0.38529827242545833, 0.019665338859497696), (0.47357895083010587, 0.031999344042706206), (0.35935583051548764, 0.0035951730371498853), (0.418774417527803, 0.021991532090154482), (0.33655979350098914, 0.013728473271936882), (0.5795742763531372, 0.006700866889811913), (0.633888494763439, 0.036280314067395784), (0.7054085401292735, 0.015219188144950999), (0.7283219380194356, 0.038371045899881606), (0.7279032948946208, 0.020915320054828775), (0.8140908945784836, 0.021769667918782345), (0.9326531697794475, 0.020552343320187443), (0.9972649387446532, 0.03218122500117304), (0.9299581288732129, 0.030817098325748242), (0.9069568906852431, 0.018944770171753167)]

        Cn6res = [(-0.04535629139876542, 0.021113893611373585), (0.0859852488954683, 0.02332645619289845), (0.10046459082681683, 0.02004970392703706), (0.132936078911145, 0.029107698275464233), (0.2173186393632259, 0.00869240606734223), (0.28203634123467103, 0.022155009014534842), (0.3643730464719448, 0.004330547417869567), (0.3626924134801637, 0.018990989593462577), (0.4647405159949755, 0.009079895858261126), (0.5022692729650203, 0.033354261771154796), (0.4694281510881381, 0.014153187156549679), (0.5566644040589314, 0.027962283011405485), (0.5451418024359513, 0.033523225772062584), (0.6813221621068155, 0.014908556442604132), (0.7433417090200445, 0.04024179292816863), (0.7711906286273695, 0.014836617946756057), (0.8194082073425, 0.029523235446106667), (0.8228541736488066, 0.0198543734727229), (0.8454411583713057, 0.017040105556118172), (0.909081676268063, 0.018714415026368024), (0.907172499406838, 0.006310711006587178)]

        alpha_u_Cn2 = [1.3973005177504485, 1.4900114242662634, 1.4161562287081038, 1.3862472778494603, 1.3712475033359504, 1.2822010660616556, 1.2476022107234472, 1.2622307499337568, 1.2129765442873435, 1.1549972762064071, 1.1349823121170068, 1.133962514187405, 1.1974279096920175, 1.15748541958521, 1.117944043265488, 1.1785343323920476, 1.1453405565825097, 1.1132352415982831, 1.102632547439785, 1.082502657644467, 1.1373164644140825]
        alpha_u_Cn4 = [1.4777834484135257, 1.4551264087420979, 1.4207589949494193, 1.3716097525234503, 1.3916646625516402, 1.3072236309459868, 1.2951939281106437, 1.2112390699261928, 1.2023504703701497, 1.130248181346157, 1.0715173721614377, 1.171582291965305, 1.1596850527927227, 1.0970915891251725, 1.1536775294698973, 1.116310670673257, 1.1236886156423231, 1.1357779370983998, 1.0867051511512427, 1.1117527004756618, 1.214311917422728]
        alpha_u_Cn6 = [1.4723108811252246, 1.5041477642870125, 1.4917032885592725, 1.452903089289555, 1.4079053024350594, 1.2838809156284545, 1.2616937264686248, 1.286026793720392, 1.201167754223943, 1.1713746818378825, 1.179985018426367, 1.1796708016489579, 1.05374140829622, 1.1316281497275797, 1.1489902037650463, 1.0977612630096436, 1.1029307965274928, 1.0954793387023247, 1.1320166475595055, 1.1146776358131572, 1.120967449783173]



    elif date == "30Apr23":
        kcs = np.arange(0.0,2.01,0.1)
        Cn2res = [(-0.11633306968492305, 0.010740171921583146), (0.02889834561250445, 0.012104671729228224), (0.14536657110775228, 0.009424554062040123), (0.23163537104547371, 0.012532305426477614), (0.15242515912558696, 0.011828347927171115), (0.3383650109102564, 0.013303459487575945), (0.32812378810336923, 0.015303794316518316), (0.3958276997010184, 0.012301162482190552), (0.5460453290916114, 0.014723038149688615), (0.5178272828040912, 0.013972401274183381), (0.5033549066549912, 0.015299011909463414), (0.6128242431292195, 0.015681116130027137), (0.6918466443805671, 0.01513918105497437), (0.6770631570702865, 0.011337144806755167), (0.6784676135416423, 0.01632924507241361), (0.7553366644956055, 0.01476431165829157), (0.7473072031259186, 0.015832862144440737), (0.7402415141386746, 0.015464888676972356), (0.10449255074334102, 0.041862147304928474), (0.9461768572345498, 0.015338113212859623), (1.0287804996912322, 0.01913611674836628)]

        alpha_u_Cn2 =  [1.4536502683367474, 1.5008864602876106, 1.4298242207760699, 1.3948865379129618, 1.3403124680111602, 1.2623773236398417, 1.3052043219371414, 1.2289067903679272, 1.2359757261695994, 1.1317656543653687, 1.2055757077200768, 1.1494504378281738, 1.1760822208883062, 1.1430329316557593, 1.1063317515497106, 1.1159600769695934, 1.1151404358811252, 1.094876277560001, 1.0223460968190237, 1.1322001802555726, 1.1723219911870746]



        Cn6res = [(-0.021496452100589723, 0.009319555659713357), (0.08380090604755229, 0.009891076560503732), (0.13970665082997513, 0.010489262712641395), (0.1112378833214494, 0.010217652125542256), (0.1897762456583451, 0.008848776057425482), (0.16656730805732914, 0.014563681607756614), (0.2985432125802804, 0.012087503461885897), (0.2552125058486404, 0.01409933642331965), (0.4275975332173634, 0.0094441211602855), (0.38440808254564895, 0.01041579451086998), (0.5599359539724217, 0.013663444884041885), (0.5979475138383012, 0.01029268204481081), (0.5598997936041182, 0.011698713178781557), (0.5860280270654185, 0.014118861948521674), (0.6131479019989021, 0.015576374473752595), (0.6780195144190486, 0.013756037790265331), (0.7184009178718741, 0.014612735157675331), (0.7451960027193545, 0.01751645464925723), (0.7282448400302478, 0.013281414637939764), (0.9064281288303504, 0.01719146097317193), (0.9679705518462316, 0.016889538681756307)]

        alpha_u_Cn6 = [1.4714654357062324, 1.4575304430869886, 1.4067202469922706, 1.3961578240297008, 1.2924201444935788, 1.3495072401951909, 1.3330172908717124, 1.2342889372070387, 1.2467632333643517, 1.186270939204213, 1.2622193959567545, 1.1208306416293747, 1.1449681730107903, 1.1417538486145935, 1.0957334853421825, 1.1143667137000213, 1.0705347911305676, 1.0932534138488819, 1.0610275596342196, 1.138363464153026, 1.0757924546471185]

        Cn10res = [(0.06268735392759966, 0.012177334177347341), (0.15994602721350706, 0.01432364634257157), (-0.02201762070093192, 0.011733974876419556), (0.0507381948970757, 0.010422982063210778), (0.12760640477598706, 0.008486235622841334), (0.14515378084939845, 0.0083004628979234), (0.21823353784378166, 0.009968221319150738), (0.26809086992378695, 0.013158923752696542), (0.3154331050104796, 0.015518491289588637), (0.45777726039912325, 0.016600252675420517), (0.3553016012221598, 0.010820693271903653), (0.6174671664618034, 0.013705834592107514), (0.5275855664364433, 0.014305117647165385), (0.5663989299235854, 0.012148950215260232), (0.6383469378789263, 0.014199791927901093), (0.6940250710098543, 0.014362743210874036), (0.8254836087904688, 0.016527532526784197), (0.7914264248253348, 0.019728166066496876), (0.8105673062023333, 0.014657977015245833), (0.7935655838902916, 0.017490215056394694), (0.8431991978341025, 0.016478712636909334)]
        #Cn10res = np.delete()

        alpha_u_Cn10 =   [100, 100, 100, 1.412646452801134, 1.3641394550966177, 1.343791066633191, 1.2346498382504412, 1.2757218181447256, 1.1515589840806653, 1.1396427962946054, 1.139608973026255, 1.1521907009625163, 1.1748493961842763, 1.0973241955431385, 1.1679076273294344, 1.1115510242401012, 1.1511004505454516, 1.080418521356136, 1.0755006110710026, 1.111711364873674, 1.0996561668047167]
        skips = {"Cn2":2,"Cn6":2,"Cn10":3}
        # old sin fit
        '''
        Cn2res = [(-0.1266917912297794, 0.028526144941257346), (0.03448459711010362, 0.021914611891257332), (0.14206239416067942, 0.017764021142885944), (0.22794243775152848, 0.032057432175757954), (0.1507376230559613, 0.007342720988852514), (0.3419144779124007, 0.017372991458620566), (0.32026997771143173, 0.023791505657524647), (0.38749059993344853, 0.022949144458460795), (0.5410837577031876, 0.026214282207906427), (0.5142430326951268, 0.021072181324653355), (0.49776940417296167, 0.017377034523919217), (0.6085401084765352, 0.014436498669558978), (0.7199551625493794, 0.04317562205592464), (0.6850685858717577, 0.04317093876367643), (0.6788219692688503, 0.014579482469368658), (0.7523789230604944, 0.01964163835734322), (0.7438863115931948, 0.021158990454350806), (0.7249128151697681, 0.022558553381528332), (0.08483305286422016, 0.1732957368363053), (0.9610505448774892, 0.041433492016104384), (1.0242414505969502, 0.018165144611271725)]

        alpha_u_Cn2 = [1.4536502683367474, 1.5008864602876106, 1.4298242207760699, 1.3948865379129618, 1.3403124680111602, 1.2623773236398417, 1.3052043219371414, 1.2289067903679272, 1.2359757261695994, 1.1317656543653687, 1.2055757077200768, 1.1494504378281738, 1.1760822208883062, 1.1430329316557593, 1.1063317515497106, 1.1159600769695934, 1.1151404358811252, 1.094876277560001, 1.0223460968190237, 1.1322001802555726, 1.1723219911870746]

        Cn6res = [(-0.01922937631836433, 0.01963807352899583), (0.08125109164756063, 0.013382281485955749), (0.12130483903370315, 0.034895253396119666), (0.12639319331359652, 0.008462851451805066), (0.2043066482110673, 0.0251905008878562), (0.1857073730147812, 0.02492102636366807), (0.2986785126450209, 0.011994533852154588), (0.25027504754513763, 0.011101630339116092), (0.43100642555692487, 0.017751279302732275), (0.3893497046413637, 0.015946610231412944), (0.555718005481073, 0.012377760605866935), (0.5878449085562648, 0.022549611536670504), (0.5557981903712348, 0.01961660997746607), (0.5896487272030112, 0.013473312549195553), (0.6100632889538995, 0.022016217928861423), (0.6783656456784983, 0.012951294059141121), (0.725789382184821, 0.01770455366515412), (0.7409865377932232, 0.032975443001650916), (0.7419725335284941, 0.018983540421732916), (0.9093475650376862, 0.028207926055303565), (0.9706653681279965, 0.013284244369809468)]

        alpha_u_Cn6 = [1.4714654357062324, 1.4575304430869886, 1.4067202469922706, 1.3961578240297008, 1.2924201444935788, 1.3495072401951909, 1.3330172908717124, 1.2342889372070387, 1.2467632333643517, 1.186270939204213, 1.2622193959567545, 1.1208306416293747, 1.1449681730107903, 1.1417538486145935, 1.0957334853421825, 1.1143667137000213, 1.0705347911305676, 1.0932534138488819, 1.0610275596342196, 1.138363464153026, 1.0757924546471185]

        Cn10res = [(0.5320036419942554, 0.16985677026033116), (0.7409560935086157, 0.08836030053038502), (-0.4555559374187401, 0.2161245170226461), (0.047344663804099595, 0.018522430562407353), (0.1252471388348498, 0.013226304820043448), (0.14305276311087384, 0.026718288029681694), (0.21325604437544712, 0.01325945077681793), (0.2620877440107482, 0.013314947076784133), (0.33569171691069577, 0.020484705398043302), (0.4575902439057407, 0.012276430931733936), (0.3907546063390977, 0.04296482633730491), (0.5922038895427746, 0.016169677009297785), (0.5387520536725244, 0.022684990772079118), (0.5659818661793371, 0.026805302422244544), (0.6281478094197216, 0.025035669793035724), (0.7039746814278789, 0.012660329725113953), (0.8272223033389206, 0.017560732609299273), (0.7879447563994834, 0.012784495007270401), (0.8099440155707626, 0.020722140839961364), (0.7953545915942407, 0.012351842158211846), (0.8300654466064288, 0.038448685070608196)]

        alpha_u_Cn10 = [-100, -100, -100, 1.412646452801134, 1.3641394550966177, 1.343791066633191, 1.2346498382504412, 1.2757218181447256, 1.1515589840806653, 1.1396427962946054, 1.139608973026255, 1.1521907009625163, 1.1748493961842763, 1.0973241955431385, 1.1679076273294344, 1.1115510242401012, 1.1511004505454516, 1.080418521356136, 1.0755006110710026, 1.111711364873674, 1.0996561668047167]
        skips = {"Cn2":2,"Cn6":1,"Cn10":4}
        '''

    elif date == "01Jul23":
        kcs = np.arange(0.0,3.01,0.1)
        Cn2res = [(0.004365314039345431, 0.01170245923677751), (0.08623858592213292, 0.010836874311700261), (0.10829717009725208, 0.012960202691998164), (0.2469638692519504, 0.011837186485035567), (0.20001431688026663, 0.015135989669358248), (0.2958014605670902, 0.013855486122443708), (0.3423241000714301, 0.014665826579437028), (0.4555394020009421, 0.017369847710818596), (0.5612610092395028, 0.014184869421358723), (0.434671657996286, 0.012865024249119788), (0.5730244494984154, 0.014238850017555022), (0.47030592854395364, 0.011651356400854695), (0.6421739023271869, 0.01794483665264144), (0.6904462665323474, 0.014976778800070891), (0.7619740022458623, 0.01719108907579778), (0.699506166360746, 0.01696035871123671), (0.8384250881052027, 0.01800269553729738), (0.800743361016663, 0.01743225937394851), (0.9452673624709222, 0.01859812076465709), (0.9930621546026797, 0.017230356764016286), (0.9802104444363773, 0.0171525830795844), (0.5195198840392945, 0.011205449396628276), (0.5527704179735323, 0.011791895404125), (0.6172319739491683, 0.011641131317195), (0.5959153112732759, 0.01178754575447469), (0.5905641965346844, 0.011905418154057289), (0.6291567856620647, 0.012435210019907203), (0.6030438276129336, 0.012637027842328173), (0.6306592532105327, 0.01146418120570706), (0.6665211538355392, 0.011610002638613604), (0.650804530873579, 0.013873967986425276)]

        alpha_u_Cn2 = [1.3975877305792235, 1.4891257509270195, 1.426630490445791, 1.3926895708658822, 1.406552341561186, 1.2780193204031083, 1.3197114226083266, 1.2672283618144509, 1.171122894404881, 1.2099100709737722, 1.1717335339241528, 1.0919199589648183, 1.1414852698937232, 1.1442691242567586, 1.1598908768369935, 1.0912569745854308, 1.1807976253204597, 1.1178056453790752, 1.2009980342212858, 1.110677414118459, 1.1365077558196097, 1.0252277497462978, 1.0257968091020797, 1.050051147462265, 0.9951514038402438, 0.9738618008953807, 1.0482199936789083, 0.9984068713298606, 1.0412087007169222, 0.9911998867948241, 0.9789783745929507]

        Cn6res = [(0.09013572895378932, 0.010703241979364555), (0.06292154672985979, 0.010258528429504655), (-0.04611559460560888, 0.009094961081250588), (0.1370248822300015, 0.008239090661638611), (0.18260113237954662, 0.010681109211578011), (0.2551471963811916, 0.008906889711364043), (0.29054816683917517, 0.01040739421508442), (0.3479840439326993, 0.012120044786998886), (0.44378729884031487, 0.013502562145026433), (0.5112319283232606, 0.012144701180684716), (0.5442606216180097, 0.013023913038323754), (0.6089224917579045, 0.010683350806377246), (0.6022458350033684, 0.017905982629205886), (0.6712586713209616, 0.015173081718932672), (0.7183009116344186, 0.016268996478780683), (0.7087793107795936, 0.014689245995774789), (0.7029917287953189, 0.01245721338743868), (0.823101916816137, 0.015260953194964486), (0.9061871314933808, 0.017090840766594846), (0.8734764372065607, 0.01759641847113761), (0.9619528566603605, 0.017007733361183743), (0.9107767527381528, 0.015601216423985665), (1.0565854099935, 0.019969422545407917), (0.6188037625021857, 0.011753816633022402), (0.5838084944917501, 0.013152675479054209), (0.6478834380006817, 0.01237072356895661), (0.5932185901140056, 0.011526981589933547), (0.6905295680752601, 0.010911198361494777), (0.6429991532590857, 0.014233619105753003), (0.6362400476882664, 0.014750402772876207), (0.6285118899368676, 0.012649111402333083)]

        alpha_u_Cn6 = [1.4556559720022126, 1.482306718945567, 1.4298160117829868, 1.3474218060873235, 1.3232965471820493, 1.3192373001797846, 1.2879301825730685, 1.265717536221737, 1.210879602720857, 1.2063260542953114, 1.1728384474361646, 1.1616309972648673, 1.135509418841953, 1.1264296894090189, 1.089933405802011, 1.1616807862760783, 1.0722849723294579, 1.1444093036207328, 1.1638791915298365, 1.1418921578600125, 1.1719769871129937, 1.153096668265857, 1.1441720437531704, 0.9985957421033524, 0.9818972689912074, 1.0380787391691766, 1.0141578294818403, 1.0398717599389362, 1.0086456574196248, 1.0059136306192562, 0.9708735665376236]

        Cn10res = [(-0.027001953371452725, 0.011933590328785815), (-0.24527447586879716, 0.01334800682153187), (0.04357825993515264, 0.015441854613581423), (0.05526131233662919, 0.009039826462550513), (0.16733055407303676, 0.010137067446978742), (0.19151638621684874, 0.013509835140241424), (0.18106627339517645, 0.009092265450395287), (0.34597624518772185, 0.010754881840447731), (0.2724531356722729, 0.0096205227645062), (0.3372786190568227, 0.011364915293406547), (0.3156500253359228, 0.012607080254048535), (0.5764290961936087, 0.017204537521425596), (0.5409569068682347, 0.010615681338837625), (0.5410701082008865, 0.014712594938754304), (0.6524027530106963, 0.014791823115756365), (0.6507141082761699, 0.014102861574435934), (0.7856399040389592, 0.01363357910408554), (0.6739785510472344, 0.01682413950137327), (0.7098870995333068, 0.015969534745753383), (0.8785907632865271, 0.018985909751451473), (0.7751074766674149, 0.015251978936365044), (0.8880176131067687, 0.015515167337124434), (0.907684729639926, 0.02103207723493013), (0.522862750620352, 0.02699820238663999), (0.5583726420853888, 0.00934184998965806), (0.5975316059860281, 0.013358654971095502), (0.5951789787512559, 0.01147484232727264), (0.6236628868838342, 0.01178368518971175), (0.6301395576682137, 0.012499002561445467), (0.6369947449000074, 0.01246681715360543), (0.6060857995549638, 0.012824693845182844)]

        alpha_u_Cn10 =  [0, 0, 1.4146812004170597, 1.3984446848237055, 1.3790721730205158, 1.2477881360443714, 1.2641193827530497, 1.2417756963325677, 1.183022427317799, 1.137037233180527, 1.148312023090561, 1.155922044602748, 1.1642050288359012, 1.1678483490035825, 1.0907488243658319, 1.080576367388524, 1.095428747826568, 1.0777070219195695, 1.0891446573764307, 1.1741190568255335, 1.0657645983798454, 1.1231575968948126, 1.1369092815732935, 1.1107554316356885, 1.012182296716012, 1.033179084602568, 0.9766338621115502, 1.0077416355484272, 0.9944553634180706, 0.9992292694651812, 0.9918021154088847]

        skips = {"Cn2":1,"Cn6":1,"Cn10":2}
        # need to check skips
    return [kcs, Cn2res,Cn6res,Cn10res,alpha_u_Cn2,alpha_u_Cn6,alpha_u_Cn10,skips]


# no longer used
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
    # remove -3th for Cn2 due to three pi walls
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


def wall_pitch_q_analysis_plot(LineWidth, FontSize, LabelSize):
    print("👌 only config + (multi mode) sinuous fit")
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.3))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")


    axcfg = plt.subplot2grid((5, 3), (0, 0), rowspan=2)
    axunphi0 = plt.subplot2grid((5, 3), (2, 0), rowspan=3)
    divider = make_axes_locatable(axunphi0)
    axunphi1 = divider.append_axes("top",size="100%",pad=0.03)
    axunphi2 = divider.append_axes("top",size="100%",pad=0.03)
    axunphi3 = divider.append_axes("top",size="100%",pad=0.03)
    axunphi4 = divider.append_axes("top",size="100%",pad=0.03)
    axunphis = [axunphi0,axunphi1,axunphi2,axunphi3,axunphi4]
    #for i in range(5):
        #axunphis.append()

    m3axcfg = plt.subplot2grid((5, 3), (0, 1), rowspan=2)
    m3axunphi0 = plt.subplot2grid((5, 3), (2, 1), rowspan=3,sharex=axunphi0,sharey = axunphi0)
    m3divider = make_axes_locatable(m3axunphi0)
    m3axunphi1 = m3divider.append_axes("top",size="100%",pad=0.03)
    m3axunphi2 = m3divider.append_axes("top",size="100%",pad=0.03)
    m3axunphi3 = m3divider.append_axes("top",size="100%",pad=0.03)
    m3axunphi4 = m3divider.append_axes("top",size="100%",pad=0.03)
    m3axunphis = [m3axunphi0,m3axunphi1,m3axunphi2,m3axunphi3,m3axunphi4]

    Caxcfg = plt.subplot2grid((5, 3), (0, 2), rowspan=2)
    Caxunphi0 = plt.subplot2grid((5, 3), (2, 2), rowspan=3,sharex=axunphi0,sharey = axunphi0)
    Cdivider = make_axes_locatable(Caxunphi0)
    Caxunphi1 = Cdivider.append_axes("top",size="100%",pad=0.03)
    Caxunphi2 = Cdivider.append_axes("top",size="100%",pad=0.03)
    Caxunphi3 = Cdivider.append_axes("top",size="100%",pad=0.03)
    Caxunphi4 = Cdivider.append_axes("top",size="100%",pad=0.03)
    Caxunphis = [Caxunphi0,Caxunphi1,Caxunphi2,Caxunphi3,Caxunphi4]

    msize = 4

    # single tube plot
    #fname = "../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn4.0.csv"
    #fname = "../data/Ne2/Dec7_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn4.0_id0.csv"
    #fname = "../data/Ne2/Dec7_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.5_Cn6.0_id0.csv"

    '''
    fname = "../data/Ne2/Apr30_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.6_Cn2.0_id0.csv"
    Cfname = "../data/Ne2/Apr30_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.6_Cn10.0_id0.csv"
    m3fname = "../data/Ne2/data_2022/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q2.5_Cn4.0.csv"
    '''

    '''
    fname = "../data/Ne2/Jun29_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.5_Cn2.0_id0.csv"
    Cfname = "../data/Ne2/Jun29_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.5_Cn10.0_id0.csv"
    m3fname = "../data/Ne2/Jun29_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q2.5_Cn2.0_id0.csv"
    '''

    fname = "../data/Ne2/Jul1_2023/data_use/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.5_Cn2.0_id1.csv"
    Cfname = "../data/Ne2/Jul1_2023/data_use/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.5_Cn10.0_id1.csv"
    m3fname = "../data/Ne2/Jul1_2023/data_use/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q2.5_Cn2.0_id1.csv"


    # note on fitting parameters:
    # C=2.0
    # [lambda_phi, phi_0] = [ 3.70942418 -1.07228748], [ 1.23443782 -0.86462099],[ 2.02118749 -0.6631666 ],  [ 1.87194865 -0.29302407],  [ 6.02565435 -0.07906099]
    # C=10.0
    # [lamda_phi, phi_0] =  [ 0.62228951 -0.85362876], [ 0.43715192 -0.6679593 ], [ 0.62382063 -0.50554623], [ 0.36347581 -0.28858196], [ 0.48607145 -0.1689633 ]

    # Jun29 data fit param

    #lam_phiC2 = np.array([ 2.40282329,2.38010938,6.31473441,2.49960228,2.12647081])
    #lam_phiC10 = np.array([0.5863347,0.54780972,0.49745928,0.50365376,0.7126027])
    #lam_phim3 = np.array([3.27391471,7.17512715,4.89385613])

    # Jul1 data fit param
    lam_phiC2 = np.array([10.50793924,])
    lam_phiC10 = np.array([])
    lam_phim3 = np.array([])


    print("C2 lam_phi:", lam_phiC2)
    print("C2 np.average(lam_phi), np.std(lam_phi):", np.average(lam_phiC2), np.std(lam_phiC2))
    print("C2 Exp(1/lam_phi):", np.exp(1/lam_phiC2))
    print("C2 Exp(-1/lam_phi):", np.exp(-1/lam_phiC2))
    print("C2 Exp(1/lam_phi)- Exp(-1/lam_phi):", np.exp(1/lam_phiC2) - np.exp(-1/lam_phiC2))
    print("-----------------------")
    print("C10 lam_phi:", lam_phiC10)
    print("C10 np.average(lam_phi), np.std(lam_phi):", np.average(lam_phiC10), np.std(lam_phiC10))
    print("C10 Exp(1/lam_phi):", np.exp(1/lam_phiC10))
    print("C10 Exp(1/lam_phi) - Exp(-1/lam_phi):", np.exp(1/lam_phiC10)- np.exp(-1/lam_phiC10))
    print("-----------------------")
    print("m3 lam_phi:", lam_phim3)
    print("m3 np.average(lam_phi), np.std(lam_phi):", np.average(lam_phim3), np.std(lam_phim3))
    print("m3 Exp(1/lam_phi):", np.exp(1/lam_phim3))
    print("m3 Exp(1/lam_phi) - Exp(-1/lam_phi):", np.exp(1/lam_phim3)- np.exp(-1/lam_phim3))
    print("-----------------------")

    ax_pitch_phi_z_plot(fname, axcfg, axunphis, 0, msize, LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4))
    axcfg.margins(y=0)
    # config plot
    x1, y1 = 0.8, 1.2+0.075/6*5
    axcfg.text(x1, y1, r"(a)", fontsize=FontSize, transform=axunphis[4].transAxes)

    ax_pitch_phi_z_plot(m3fname, m3axcfg, m3axunphis, 0, msize, LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4), m=3)
    m3axcfg.margins(y=0)
    # config plot
    x1, y1 = 0.8, 1.2+0.075/6*5
    m3axcfg.text(x1, y1, r"(b)", fontsize=FontSize, transform=m3axunphis[4].transAxes)

    ax_pitch_phi_z_plot(Cfname, Caxcfg, Caxunphis, 0, msize, LineWidth, nbin=5, z_relative_range=(1 / 4, 3 / 4))
    Caxcfg.margins(y=0)
    # config plot
    x1, y1 = 0.8, 1.2+0.075/6*5
    Caxcfg.text(x1, y1, r"(c)", fontsize=FontSize, transform=Caxunphis[4].transAxes)





    # print(filename, params[1])
    # axs[2].scatter(phi,x)
    for i in range(5):
        axunphis[i].tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=False, labelleft=True, labelsize=LabelSize)
        #axunphis[i].set_xlim(-2.7, 2.7)
        axunphis[i].yaxis.set_major_locator(MultipleLocator(0.5))
        axunphis[i].yaxis.set_minor_locator(MultipleLocator(0.25))
        axunphis[i].set_ylim(-0.1, 1.1)
    #axunphi0.set_xlabel(r"$\phi = \arctan(y/x)$",fontsize=FontSize)
    axunphi0.xaxis.set_major_locator(MultipleLocator(2))
    axunphi0.xaxis.set_minor_locator(MultipleLocator(1))
    axunphi0.tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=True, labelleft=True, labelsize=LabelSize)
    #axunphis[2].set_ylabel(r"$\left<(\vu{u}\cdot \vu{n})^2\right>$",fontsize=FontSize)
    axunphis[2].set_ylabel(r"$(\vu{u}\cdot \vu{n})^2$",fontsize=FontSize)

    x1, y1 = 0.09, 0.057*5/9*5
    axunphi0.text(x1, y1, r"(d)", fontsize=FontSize, transform=axunphis[0].transAxes)

    for i in range(5):
        m3axunphis[i].tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=False, labelleft=False, labelsize=LabelSize)
        #axunphis[i].set_xlim(-2.7, 2.7)
        m3axunphis[i].yaxis.set_major_locator(MultipleLocator(0.5))
        m3axunphis[i].yaxis.set_minor_locator(MultipleLocator(0.25))
        m3axunphis[i].set_ylim(-0.1, 1.1)
    m3axunphi0.set_xlabel(r"$\phi = \arctan(y/x)$",fontsize=FontSize)
    m3axunphi0.xaxis.set_major_locator(MultipleLocator(2))
    m3axunphi0.xaxis.set_minor_locator(MultipleLocator(1))
    m3axunphi0.tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=True, labelleft=False, labelsize=LabelSize)
    #m3axunphis[2].set_ylabel(r"$(\vu{u}\cdot \vu{n})^2$",fontsize=FontSize)
    x1, y1 = 0.09, 0.057*5/9*5
    m3axunphi0.text(x1, y1, r"(e)", fontsize=FontSize, transform=m3axunphis[0].transAxes)

    for i in range(5):
        Caxunphis[i].tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=False, labelleft=False, labelsize=LabelSize)
        #axunphis[i].set_xlim(-2.7, 2.7)
        Caxunphis[i].yaxis.set_major_locator(MultipleLocator(0.5))
        Caxunphis[i].yaxis.set_minor_locator(MultipleLocator(0.25))
        Caxunphis[i].set_ylim(-0.1, 1.1)
    #Caxunphi0.set_xlabel(r"$\phi = \arctan(y/x)$",fontsize=FontSize)
    Caxunphi0.xaxis.set_major_locator(MultipleLocator(2))
    Caxunphi0.xaxis.set_minor_locator(MultipleLocator(1))
    Caxunphi0.tick_params(which="both", direction="in", bottom="on", top="on", right="on", left="on", labelbottom=True, labelleft=False, labelsize=LabelSize)
    #Caxunphis[2].set_ylabel(r"$(\vu{u}\cdot \vu{n})^2$",fontsize=FontSize)
    x1, y1 = 0.09, 0.057*5/9*5
    Caxunphi0.text(x1, y1, r"(f)", fontsize=FontSize, transform=Caxunphis[0].transAxes)



    #plt.tight_layout(pad=0.1)
    plt.subplots_adjust(wspace=0.05,hspace=0.3,left=0.12,right=0.97,bottom=0.08,top=0.98)
    plt.savefig("figures/wall_pitch_analysis.pdf", format="pdf")
    plt.close()



def wall_pitch_q_result_plot(LineWidth, FontSize, LabelSize):
    print("👌")
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axtankc = plt.subplot2grid((1, 2), (0, 0))
    axthekc = plt.subplot2grid((1, 2), (0, 1), sharex = axtankc)

    msize = 4
    # no longer needed, since R already considered in fitting slope
    #kcs, Cn2res,Cn6res,Cn10res,alpha_u_Cn2,alpha_u_Cn6,alpha_u_Cn10, skips = wall_pitch_tan_data_get(date = "30Apr23")
    kcs, Cn2res,Cn6res,Cn10res,alpha_u_Cn2,alpha_u_Cn6,alpha_u_Cn10, skips = wall_pitch_tan_data_get(date = "01Jul23")
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
    #axtankc.errorbar(np.delete(kcs[n0Cn2:],-3),np.delete(Cn2ta[n0Cn2:],-3),yerr = np.delete(Cn2taerr[n0Cn2:],-3), color = "b", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=2")
    axtankc.errorbar(kcs[n0Cn2:],Cn2ta[n0Cn2:],yerr = Cn2taerr[n0Cn2:], color = "b", marker = "o", linestyle="None" ,ms=msize,mfc="None", label="C=2")
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
    axtankc.set_ylim(-0.01,1.151)
    axtankc.legend(loc="upper left", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)

    x1, y1 = 0.75, 0.075
    axtankc.text(x1, y1, r"(a)", fontsize=FontSize, transform=axtankc.transAxes)


    axthekc.plot(kcs[n0Cn2:],alpha_u_Cn2[n0Cn2:] - np.arctan(Cn2ta[n0Cn2:]), color = "b", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=2")
    axthekc.plot(kcs[n0Cn6:],alpha_u_Cn6[n0Cn6:] - np.arctan(Cn6ta[n0Cn6:]), color = "r", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=6")
    axthekc.plot(kcs[n0Cn10:],alpha_u_Cn10[n0Cn10:] - np.arctan(Cn10ta[n0Cn10:]), color = "k", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=10")
    axthekc.set_xlabel(r"$k_c$",fontsize=FontSize)
    axthekc.set_ylabel(r"$\theta = \alpha_u - \alpha$",fontsize=FontSize)

    # following are for plotting aphla_u, just for checking
    '''
    axthekc.plot(np.delete(kcs[n0Cn2:],-3),np.delete(alpha_u_Cn2[n0Cn2:],-3) , color = "b", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=2")
    axthekc.plot(kcs[n0Cn6:],alpha_u_Cn6[n0Cn6:], color = "r", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=6")
    axthekc.plot(kcs[n0Cn10:],alpha_u_Cn10[n0Cn10:], color = "k", marker = "s", linestyle="None", ms=msize, mfc = "None", label = "C=10")
    axthekc.set_xlabel(r"$k_c$",fontsize=FontSize)
    axthekc.set_ylabel(r"$\alpha_u$",fontsize=FontSize)
    '''
    axthekc.yaxis.set_label_position("left")
    #axthekc.yaxis.tick_right()
    axthekc.tick_params(which="both", direction="in", bottom="on", top="off", right="off", left="on", labelbottom=True, labelleft=True, labelsize=LabelSize)

    axthekc.xaxis.set_major_locator(MultipleLocator(1.0))
    axthekc.xaxis.set_minor_locator(MultipleLocator(0.2))
    axthekc.yaxis.set_major_locator(MultipleLocator(0.5))
    axthekc.yaxis.set_minor_locator(MultipleLocator(0.1))
    axthekc.set_ylim(-0.01,1.51)
    #axthekc.set_ylim(1,1.51)
    axthekc.legend(loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    axthekc.text(0.75, 0.075, r"(b)", fontsize=FontSize, transform=axthekc.transAxes)

    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.4,hspace=0.3,left=0.12,right=0.97,bottom=0.2,top=0.98)
    plt.savefig("figures/wall_pitch_result.pdf", format="pdf")
    plt.close()