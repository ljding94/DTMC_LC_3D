import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.colors import Normalize

def n(m, alpha, gamma, phi, z):
    nx = np.cos(phi)
    ny = np.sin(phi)
    nz = np.zeros(len(phi))
    return np.array([nx, ny, nz])


def u(m, alpha, gamma, phi, z):
    ux = (np.cos(phi) * np.cos((m * (phi - z * np.tan(alpha))) / 2.0) - (-1 + gamma + gamma * np.sin(alpha)) * np.sin(phi) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    uy = (np.cos((m * (phi - z * np.tan(alpha))) / 2.0) * np.sin(phi) + np.cos(phi) * (-1 + gamma + gamma * np.sin(alpha)) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    uz = (gamma * np.cos(alpha) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    return np.array([ux, uy, uz])


def beta_cal(m,alpha,phi,z):
    return m/2*(phi-np.tan(alpha)*z)


def u_beta(m, alpha, gamma, phi, z):
    beta = beta_cal(m,alpha,phi,z)
    ux = (np.cos(beta)*np.cos(phi) - (-1 + gamma + gamma*np.sin(alpha))*np.sin(beta)*np.sin(phi))/np.sqrt(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)
    uy = (np.cos(phi)*(-1 + gamma + gamma*np.sin(alpha))*np.sin(beta) + np.cos(beta)*np.sin(phi))/np.sqrt(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)
    uz = (gamma*np.cos(alpha)*np.sin(beta))/np.sqrt(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)

    return np.array([ux,uy,uz])

def test():
    bn_phi,bn_z = 100,100
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 - 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 - 0.5 / bn_z, bn_z))
    phi, z = phi.flatten(), z.flatten()
    m = 2
    alpha = 0.1
    gamma = 0.1
    print(np.max(u_beta(m, alpha, gamma, phi, z) - u(m, alpha, gamma, phi, z)))
    # tested, they are the same


## sample model config plot


def ax_2mod_u_plot(ax,m,alpha,gamma,bn_phi,bn_z,d=0.1):
    bns_phi, bns_z = 100, 50
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bns_phi, 1 + 0.5 / bns_phi, bns_phi), np.linspace(0.5 / bns_z, 1 + 0.5 / bns_z, bns_z))
    x, y = np.cos(phi), np.sin(phi)

    print("np.shape(x)",np.shape(x))
    print("np.shape(y)",np.shape(y))
    print("np.shape(z)",np.shape(z))
    # plot surface mesh
    ax.plot_surface(x,y,z, linewidth=0.5,shade=0,color = "gray", edgecolor="black",alpha=0.5,rstride=250, cstride=20)
    #ax.plot_surface(x,y,z, linewidth=0.5,shade=0,color = "gray",alpha=0.5,rstride=10, cstride=10)


    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 + 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 + 0.5 / bn_z, bn_z))
    x, y = np.cos(phi), np.sin(phi)
    x,y,z,phi = x.flatten(),y.flatten(),z.flatten(),phi.flatten()
    ux,uy,uz = u(m, alpha, gamma, phi, z)
    nx,ny,nz = n(m, alpha, gamma, phi, z)
    deg=np.arccos(np.abs(ux*nx+uy*ny+uz*nz))

    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    # plot director
    for i in range(len(x)):
        #ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(abs_un[i])),label=r"$u$")
        ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=2,color=cmap(norm(deg[i])))
    #ax.set_aspect("auto")
    ax.view_init(elev=50) #., azim=-40)
    ax.dist=7.2
    ax.set_axis_off()
    ax.set_frame_on(False)



def demo_config_2mod_u_plot(LineWidth, FontSize, LabelSize):
    print("plotting sample configure demostrating combination of 2 mod of director field on cylinder")


    ppi = 72
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    fig, axs = plt.subplots(3, 3,figsize=(246 / ppi * 1, 246 / ppi * 1),subplot_kw=dict(projection='3d'))
    # Remove vertical space between axes
    #fig.subplots_adjust(hspace=-0.5)
    # Plot each graph, and manually set the y tick values
    alphas = [0,0.4,0.8]
    gammas = [0,1,0.5]
    bn_phi,bn_z = 20,10
    for i in range(3):
        axs[i,0].text2D(-0.05,0.4, r"$\alpha=%.1f$"%alphas[i], fontsize=FontSize, transform=axs[i,0].transAxes,rotation=90)
        axs[0,i].text2D(0.4,1.05, r"$\gamma=%.1f$"%gammas[i], fontsize=FontSize, transform=axs[0,i].transAxes)
        for j in range(3):
            ax_2mod_u_plot(axs[i,j],2,alphas[i],gammas[j],bn_phi,bn_z,d=0.2)
    plt.tight_layout(pad=0.37)
    #plt.tight_layout(pad=-1)
    #plt.show()
    plt.savefig("figures/twomode_config_demo.pdf", format="pdf")





def get_energy_color_map_data():
    pass

def energy_color_map_plot(LineWidth, FontSize, LabelSize):
    pass
    # can use diverging color maps
    # e.g. 'coolwarm'

