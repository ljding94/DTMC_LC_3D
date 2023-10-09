import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.colors as colors
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from scipy.signal import savgol_filter

def n(m, alpha, gamma, phi, z):
    nx = np.cos(phi)
    ny = np.sin(phi)
    nz = np.zeros(len(phi))
    return np.array([nx, ny, nz])


def u_misaligned(m, alpha, gamma, phi, z):
    ux = (np.cos(phi) * np.cos((m * (phi - z * np.tan(alpha))) / 2.0) - (-1 + gamma + gamma * np.sin(alpha)) * np.sin(phi) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    uy = (np.cos((m * (phi - z * np.tan(alpha))) / 2.0) * np.sin(phi) + np.cos(phi) * (-1 + gamma + gamma * np.sin(alpha)) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    uz = (gamma * np.cos(alpha) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    return np.array([ux, uy, uz])

def u(m, alpha, gamma, phi, z):
    ux = (np.cos(phi)*np.cos((m*(phi - z*np.tan(alpha)))/2.) + (1 - gamma + gamma*np.sin(alpha))*np.sin(phi)*np.sin((m*(phi - z*np.tan(alpha)))/2.))/np.sqrt(1 + (-1 + gamma)*gamma + (-1 + gamma)*gamma*(np.cos(m*(phi - z*np.tan(alpha)))*(-1 + np.sin(alpha)) - np.sin(alpha)))
    uy= (np.cos((m*(phi - z*np.tan(alpha)))/2.)*np.sin(phi) + np.cos(phi)*(-1 + gamma - gamma*np.sin(alpha))*np.sin((m*(phi - z*np.tan(alpha)))/2.))/np.sqrt(1 + (-1 + gamma)*gamma + (-1 + gamma)*gamma*(np.cos(m*(phi - z*np.tan(alpha)))*(-1 + np.sin(alpha)) - np.sin(alpha)))
    uz= -((gamma*np.cos(alpha)*np.sin((m*(phi - z*np.tan(alpha)))/2.))/np.sqrt(1 + (-1 + gamma)*gamma + (-1 + gamma)*gamma*(np.cos(m*(phi - z*np.tan(alpha)))*(-1 + np.sin(alpha)) - np.sin(alpha))))
    return np.array([ux,uy,uz])

def uz_at_wall(alpha,gamma):
    ans =  (gamma*np.cos(alpha))/np.sqrt(1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))
    return ans

def uz_align_at_wall(alpha,gamma):
    ans = (gamma*np.cos(alpha))/np.sqrt(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(-1 + 2*np.sin(alpha)))
    return ans

def alpha_u_at_wall(alpha,gamma):
    # use sin, whith has sign -issue, e.g. 3/4\pi ->-1/4\pi
    #ans = -(-1 + gamma + gamma*np.sin(alpha))/np.sqrt(1 - 2*gamma + 2*gamma**2 + 2*(-1 + gamma)*gamma*np.sin(alpha))
    #return np.arcsin(ans)
    #ans = -((-1 + gamma)/np.cos(alpha))/gamma + np.tan(alpha)
    ans = (gamma*np.cos(alpha))/np.sqrt(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(-1 + 2*np.sin(alpha)))
    return np.arccos(ans)



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


def ax_2mod_u_3Dplot(ax,m,alpha,gamma,d=0.1):
    bns_phi, bns_z = 101, 51
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.0 / bns_phi, 1 + 0.0 / bns_phi, bns_phi), np.linspace(0.5 / bns_z, 1 + 0.5 / bns_z, bns_z))
    x, y = np.cos(phi), np.sin(phi)

    print("np.shape(x)",np.shape(x))
    print("np.shape(y)",np.shape(y))
    print("np.shape(z)",np.shape(z))
    # plot surface mesh
    ax.plot_surface(x,y,z, linewidth=0,shade=0,color = "gray", edgecolor="gray",alpha=0.6) # ,rstride=1, cstride=1)
    #ax.plot_surface(x,y,z, linewidth=0.5,shade=0,color = "gray",alpha=0.5,rstride=10, cstride=10)

    bn_phi,bn_z = 18,9
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.0 / bn_phi, 1 + 0.0 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 + 0.5 / bn_z, bn_z))
    x, y = np.cos(phi), np.sin(phi)
    x,y,z,phi = x.flatten(),y.flatten(),z.flatten(),phi.flatten()
    ux,uy,uz = u(m, alpha, gamma, phi, z)
    #ux,uy,uz = u(m, alpha, gamma, phi, z)
    nx,ny,nz = n(m, alpha, gamma, phi, z)
    deg=np.arccos(np.abs(ux*nx+uy*ny+uz*nz))

    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    # plot director
    for i in range(len(x)):
        #ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(abs_un[i])),label=r"$u$")
        ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=2,color=cmap(norm(deg[i])))
        # for the propurse of seeing direction
    #ax.scatter([x-0.5*d*ux],[y-0.5*d*uy],[z-0.5*d*uz],s=5,marker="o",color=cmap(norm(deg)))

    #plot wall line
    z = np.linspace(0,1,50)
    delphi = np.pi/m
    while delphi<2*np.pi:
        phi = np.tan(alpha)*z + delphi
        x = np.cos(phi)
        y = np.sin(phi)
        ax.plot3D(x,y,z, linestyle="-", color="k", alpha=0.9, linewidth = 1)
        delphi += 2*np.pi/m

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
    alphas = [0,np.pi/8,np.pi/4]
    alphas_label = ["0","\pi/8","\pi/4"]
    alphas_u = []
    uzs = []
    thetas_u = []
    thetas_label = [["\pi/2","\pi/4","0"],["3\pi/8","3\pi/16","0",],["\pi/4","\pi/8","0"]]
    gammas = [0,0.5,1]
    for i in range(3):
        axs[i,0].text2D(-0.05,0.4, r"$\alpha=%s$"%alphas_label[i], fontsize=FontSize, transform=axs[i,0].transAxes,rotation=90)
        axs[0,i].text2D(0.4,1.05, r"$\gamma=%.1f$"%gammas[i], fontsize=FontSize, transform=axs[0,i].transAxes)
        for j in range(3):
            ax_2mod_u_3Dplot(axs[i,j],3,alphas[i],gammas[j],d=0.2)
            #uz = uz_at_wall(alphas[i],gammas[j])
            alpha_u = alpha_u_at_wall(alphas[i],gammas[j])
            #if(alpha_u<0):
            #    alpha_u = np.pi+alpha_u
            alphas_u.append(alpha_u)
            thetas_u.append(alpha_u-alphas[i])
            print(alpha_u/np.pi*16)
            #axs[i,j].text2D(0.3,-0.05,r"$\theta=%s$"%thetas_label[i][j],fontsize=FontSize, transform=axs[i,j].transAxes)
    print("alphas_u/np.pi*16",np.array(alphas_u)/np.pi*16)
    print("(thetas_u)/np.pi*16",np.array(thetas_u)/np.pi*16)
    plt.tight_layout(pad=0.37)
    #plt.tight_layout(pad=-1)
    #plt.show()
    plt.savefig("figures/twomode_config_demo_m3.pdf", format="pdf")

def ax_2mod_2Dplot(ax,mod,view,hshift=0,vshift=0,du=0.2,dt=0.1,label="",FontSize=9):
    # get ux,uy,uz for different mod
    bn_phi,bn_z = 16,1
    #phi, z = np.meshgrid(2 * np.pi * np.linspace(0.0 / bn_phi, 1 + 0.0 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 + 0.5 / bn_z, bn_z))
    phi, z = np.arange(0,2*np.pi,2*np.pi/bn_phi)+2*np.pi/bn_phi,np.zeros(bn_phi)
    print(np.shape(phi))
    x, y = np.cos(phi), np.sin(phi)
    x,y,z,phi = x.flatten(),y.flatten(),z.flatten(),phi.flatten()

    nx,ny,nz = n(2, 0, 0, phi, z)

    if(mod=="nematic-z"):
        ux,uy,uz = np.zeros(len(x)),np.zeros(len(x)),np.ones(len(x))
    elif(mod=="nematic-x"):
        ux,uy,uz = u(2, 0, 0, phi, z)
    elif(mod=="smectic"):
        ux,uy,uz = u(0, 0, 0, phi, z)
    elif(mod=="cholesteric"):
        ux,uy,uz = u(2, 0, 1, phi, z)
    else:
        print("\nundefined mod\n")


    print("ux:",ux)
    print("uy:",uy)
    deg=np.arccos(np.abs(ux*nx+uy*ny+uz*nz))
    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)
    print("np.min(deg),np.max(deg)",np.min(norm(deg)),np.max(norm(deg)))
    vlw = 1.5
    ms = 3
    if(view=="xy"):
        x = x+ hshift
        y = y + vshift
        #plot circle
        pphi = np.linspace(0,2*np.pi,200)
        ax.plot(np.cos(pphi)+hshift,np.sin(pphi)+vshift,"-",color="gray",lw=1)
        ax.text(-1.4+hshift,0.7+vshift,label,fontsize=FontSize)
        #plot director in T shape
        print(mod,"xy len(x)",len(x),x,ux)
        for i in range(len(x)):
            print("i",i)
            ax.plot([x[i],x[i]+du*ux[i]],[y[i],y[i]+du*uy[i]],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            #ax.plot([x[i]-0.5*du*ux[i],x[i]+0.5*du*ux[i]],[y[i]-0.5*du*uy[i],y[i]+0.5*du*uy[i]],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            #ax.text(x[i],y[i],i)
            if(abs(uz[i])>1-1e-4):
                ax.scatter([x[i]],[y[i]],color=cmap(norm(deg[i])),s=ms,zorder=10)
            if(abs(uz[i])>1e-4 and abs(uz[i])<(1-1e-4)):
                tx,ty,tz = np.cross([ux[i],uy[i],uz[i]],[0,0,1])
                mag_t = np.sqrt(np.inner([tx,ty,tz],[tx,ty,tz]))
                if(mag_t==0):
                    tx,ty,tz = np.sin(phi[i]),-np.cos(phi[i]),0
                else:
                    tx,ty,tz=tx/mag_t,ty/mag_t,tz/mag_t
                print("inner t=", np.inner([tx,ty,tz],[tx,ty,tz]))
                print("mag_t=",mag_t)
                # potition of tilt is the end of the rods most outwards
                direct = np.sign(uz[i])/2+0.5
                direct *= du
                ax.plot([x[i]+direct*ux[i]-0.5*dt*tx,x[i]+direct*ux[i]+0.5*dt*tx],[y[i]+direct*uy[i]-0.5*dt*ty,y[i]+direct*uy[i]+0.5*dt*ty],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            #ax.plot([x[i]-0.5*du*ux[i]-0.5*dt*tx,x[i]-0.5*du*ux[i]+0.5*dt*tx],[y[i]-0.5*du*uy[i]-0.5*dt*ty,y[i]-0.5*du*uy[i]+0.5*dt*ty],"-",linewidth=vlw,color=cmap(norm(deg[i])))

    elif(view=="xz"):
        x = x+ hshift
        z = z + vshift
        #plot box
        zshift=-1.6
        ax.plot(np.array([-1,1,1,-1,-1])+hshift,0.5*(np.array([-0.5,-0.5,0.5,0.5,-0.5]))+vshift+zshift,"-",color="gray",lw=1)
        #ax.text(1.4+hshift,-0.1+vshift+zshift,label,fontsize=FontSize)
        #plot director in T shape
        for i in range(-bn_phi//2-1,0):
            #ax.text(x[i],z[i]+zshift,i)
            ax.plot([x[i],x[i]+du*ux[i]],[z[i]+zshift,z[i]+zshift+du*uz[i]],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            #ax.plot([x[i]-0.5*du*ux[i],x[i]+0.5*du*ux[i]],[z[i]+zshift-0.5*du*uz[i],z[i]+zshift+0.5*du*uz[i]],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            if(abs(uy[i])>1-1e-4):
                ax.scatter([x[i]],[z[i]+zshift],color=cmap(norm(deg[i])),s=ms,zorder=10)
            if(abs(uy[i])>1e-4 and abs(uy[i])<(1-1e-4)):
                tx,ty,tz = np.cross([ux[i],uy[i],uz[i]],[0,1,0])
                mag_t = np.sqrt(np.inner([tx,ty,tz],[tx,ty,tz]))
                tx,ty,tz=tx/mag_t,ty/mag_t,tz/mag_t
                '''
                if(mag_t==0):
                    tx,ty,tz = np.sin(phi[i]),-np.cos(phi[i]),0
                else:
                    tx,ty,tz=tx/mag_t,ty/mag_t,tz/mag_t
                '''
                direct = -np.sign(uy[i])/2+0.5
                direct *= du
                ax.plot([x[i]+direct*ux[i]-0.5*dt*tx,x[i]+direct*ux[i]+0.5*dt*tx],[z[i]+direct*uz[i]+zshift-0.5*dt*tz,z[i]+direct*uz[i]+zshift+0.5*dt*tz],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            #ax.plot([x[i]-0.5*du*ux[i]-0.5*dt*tx,x[i]-0.5*du*ux[i]+0.5*dt*tx],[z[i]+zshift-0.5*du*uz[i]-0.5*dt*tz,z[i]+zshift-0.5*du*uz[i]+0.5*dt*tz],"-",linewidth=vlw,color=cmap(norm(deg[i])))
    elif(view=="phiz"):
        uphi = -np.sin(phi)*ux+np.cos(phi)*uy
        phi = phi - 3*np.pi/2 # take phi \in (pi,2pi), center it
        phi = phi + hshift
        z = z + vshift


        nun = -(ux*nx+uy*ny+uz*nz) # inwards normal direction, y-like for phi-z
        # map (phi->x) (nun->y)
        print("np.shape(uphi),np.shape(nun)",np.shape(uphi),np.shape(nun),np.shape(z),np.shape(uz),np.shape(deg))
        print("uphi[0]",uphi[0])

        #plot box
        zshift=-2.4
        ax.plot(np.pi/2*np.array([-1,1,1,-1,-1])+hshift,0.5*(np.array([-0.5,-0.5,0.5,0.5,-0.5]))+vshift+zshift,"-",color="gray",lw=1)
        ax.text(1.1+hshift,-0.5+vshift+zshift,label,fontsize=FontSize)
        #plot director in T shape
        for i in range(-bn_phi//2-1,0):
            #ax.text(phi[i],z[i]+zshift,i)
            ax.plot([phi[i],phi[i]+du*uphi[i]],[z[i]+zshift,z[i]+zshift+du*uz[i]],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            if(abs(nun[i])>1-1e-4):
                ax.scatter([phi[i]],[z[i]+zshift],color=cmap(norm(deg[i])),s=ms,zorder=10)
            if(abs(nun[i])>1e-4 and abs(nun[i])<(1-1e-4)):
                tx,ty,tz = np.cross([uphi[i],nun[i],uz[i]],[0,1,0])
                mag_t = np.sqrt(np.inner([tx,ty,tz],[tx,ty,tz]))
                tx,ty,tz=tx/mag_t,ty/mag_t,tz/mag_t
                '''
                if(mag_t==0):
                    tx,ty,tz = np.sin(phi[i]),-np.cos(phi[i]),0
                else:
                    tx,ty,tz=tx/mag_t,ty/mag_t,tz/mag_t
                '''
                direct = -np.sign(nun[i])/2+0.5
                direct *= du
                ax.plot([phi[i]+direct*uphi[i]-0.5*dt*tx,phi[i]+direct*uphi[i]+0.5*dt*tx],[z[i]+direct*uz[i]+zshift-0.5*dt*tz,z[i]+direct*uz[i]+zshift+0.5*dt*tz],"-",linewidth=vlw,color=cmap(norm(deg[i])))
    elif(view=="phinrho"):
        uphi = -np.sin(phi)*ux+np.cos(phi)*uy
        unrho = -(ux*nx+uy*ny+uz*nz) # inwards normal direction, y-like for phi-z
        phi = phi - 3*np.pi/2 # take phi \in (pi,2pi), center it
        phi = phi + hshift
        nrho = np.zeros(len(phi)) + vshift
        # negative rho direction
        #plot box
        nrhoshift=-3.1
        #ax.plot(np.pi/2*np.array([-1,1,1,-1,-1])+hshift,0.1*(np.array([-0.5,-0.5,0.5,0.5,-0.5]))+vshift+nrhoshift,"-",color="gray",lw=1)
        ax.plot(np.pi/2*np.array([-1,1])+hshift,0.1*(np.array([0,0]))+vshift+nrhoshift,"-",color="gray",lw=1)
        #plot director in T shape
        for i in range(-bn_phi//2-1,0):
            #ax.text(phi[i],z[i]+zshift,i)
            ax.plot([phi[i],phi[i]+du*uphi[i]],[nrho[i]+nrhoshift,nrho[i]+nrhoshift+du*unrho[i]],"-",linewidth=vlw,color=cmap(norm(deg[i])))
            if(abs(uz[i])>1-1e-4):
                ax.scatter([phi[i]],[nrho[i]+nrhoshift],color=cmap(norm(deg[i])),s=ms,zorder=10)
            if(abs(uz[i])>1e-4 and abs(uz[i])<(1-1e-4)):
                tx,ty,tz = np.cross([uphi[i],unrho[i],uz[i]],[0,0,1])
                mag_t = np.sqrt(np.inner([tx,ty,tz],[tx,ty,tz]))
                tx,ty,tz=tx/mag_t,ty/mag_t,tz/mag_t
                '''
                if(mag_t==0):
                    tx,ty,tz = np.sin(phi[i]),-np.cos(phi[i]),0
                else:
                    tx,ty,tz=tx/mag_t,ty/mag_t,tz/mag_t
                '''
                direct = np.sign(uz[i])/2+0.5
                direct *= du
                ax.plot([phi[i]+direct*uphi[i]-0.5*dt*tx,phi[i]+direct*uphi[i]+0.5*dt*tx],[nrho[i]+direct*unrho[i]+nrhoshift-0.5*dt*tz,nrho[i]+direct*unrho[i]+nrhoshift+0.5*dt*tz],"-",linewidth=vlw,color=cmap(norm(deg[i])))


from mpl_toolkits.axes_grid1.inset_locator import inset_axes
def two_mod_diagram_plot(LineWidth, FontSize, LabelSize):
    # diagram for explaining monte carlo results: alignment of directors
    # plot sample config of nematic and cholesteric walls using Kleman's nail representation
    print("plotting illustration diagrams for walls")

    ppi = 72
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    fig, ax = plt.subplots(1, 1,figsize=(246 / ppi * 1, 246 / ppi * 1.05) )

    delh,delv=3.6,4.6
    delh0 = 1.0

    ax_2mod_2Dplot(ax,"smectic","xy",hshift=0,vshift=delv,label=r"(a)")
    ax_2mod_2Dplot(ax,"smectic","xz",hshift=0,vshift=delv)
    ax_2mod_2Dplot(ax,"smectic","phiz",hshift=0,vshift=delv)
    ax_2mod_2Dplot(ax,"smectic","phinrho",hshift=0,vshift=delv)

    ax_2mod_2Dplot(ax,"nematic-x","xy",hshift=delh,vshift=delv,label=r"(b)")
    ax_2mod_2Dplot(ax,"nematic-x","xz",hshift=delh,vshift=delv)
    ax_2mod_2Dplot(ax,"nematic-x","phiz",hshift=delh,vshift=delv)
    ax_2mod_2Dplot(ax,"nematic-x","phinrho",hshift=delh,vshift=delv)

    ax_2mod_2Dplot(ax,"cholesteric","xy",hshift=0,vshift=0,label=r"(c)")
    ax_2mod_2Dplot(ax,"cholesteric","xz",hshift=0,vshift=0)
    ax_2mod_2Dplot(ax,"cholesteric","phiz",hshift=0,vshift=0)
    ax_2mod_2Dplot(ax,"cholesteric","phinrho",hshift=0,vshift=0)

    ax_2mod_2Dplot(ax,"nematic-z","xy",hshift=delh,vshift=0,label=r"(d)")
    ax_2mod_2Dplot(ax,"nematic-z","xz",hshift=delh,vshift=0)
    ax_2mod_2Dplot(ax,"nematic-z","phiz",hshift=delh,vshift=0)
    ax_2mod_2Dplot(ax,"nematic-z","phinrho",hshift=delh,vshift=0)


    axins = ax.inset_axes([0.02,0.05,0.02,0.35])

    cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),cax=axins,ticks=[0,np.pi/6,np.pi/3,np.pi/2])
    cbar.ax.set_yticklabels([r"$0$",r"$\pi/6$",r"$\pi/3$",r"$\pi/2$"],fontsize=FontSize)
    cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    cbar.ax.set_title(r"$\arccos{|\vu{u}\cdot\vu{n}|}$",fontsize=FontSize)


    #axcoor = ax.inset_axes([-0.1,0.55,0.02,0.35])

    # plot the follorgin coordinate in axcoor
    delar = 0.4
    delh_mt = -0.7
    delvar = delv
    lshift = - 0.27
    hw = 0.05
    '''
    axcoor.arrow(delh_mt*delh+0.5,delvar, delar,0, head_width=hw)
    axcoor.arrow(delh_mt*delh+0.5,delvar, 0,delar, head_width=hw)
    axcoor.text(delh_mt*delh+delar+0.5,delvar+lshift, r"$x$", fontsize=FontSize)
    axcoor.text(delh_mt*delh+lshift+0.5,delvar+delar, r"$y$", fontsize=FontSize)

    axcoor.arrow(delh_mt*delh+0.3,-1.7+delvar, delar,0, head_width=hw)
    axcoor.arrow(delh_mt*delh+0.3,-1.7+delvar, 0,delar, head_width=hw)
    axcoor.text(delh_mt*delh+delar+0.3,+lshift-1.7+delvar, r"$x$", fontsize=FontSize)
    axcoor.text(delh_mt*delh+lshift+0.3,delar-1.7+delvar, r"$z$", fontsize=FontSize)

    axcoor.arrow(delh_mt*delh,-2.55+delvar, delar,0, head_width=hw)
    axcoor.arrow(delh_mt*delh,-2.55+delvar, 0,delar, head_width=hw)
    axcoor.text(delh_mt*delh+delar-lshift,lshift-2.55+delvar, r"$\phi$", fontsize=FontSize)
    axcoor.text(delh_mt*delh+lshift,delar-2.55+delvar, r"$z$", fontsize=FontSize)

    axcoor.arrow(delh_mt*delh,-2.75+delvar, delar,0, head_width=hw)
    axcoor.arrow(delh_mt*delh,-2.75+delvar, 0,-delar, head_width=hw)
    #axcoor.text(delh_mt*delh+delar,+lshift-2.95+delvar, r"$\phi$", fontsize=FontSize)
    axcoor.text(delh_mt*delh+lshift,-delar-2.75+delvar, r"$\rho$", fontsize=FontSize)

    '''
    ax.arrow(delh_mt*delh+0.5,delvar, delar,0, head_width=hw)
    ax.arrow(delh_mt*delh+0.5,delvar, 0,delar, head_width=hw)
    ax.text(delh_mt*delh+delar+0.5,delvar+lshift, r"$x$", fontsize=FontSize)
    ax.text(delh_mt*delh+lshift+0.5,delvar+delar, r"$y$", fontsize=FontSize)

    ax.arrow(delh_mt*delh+0.3,-1.7+delvar, delar,0, head_width=hw)
    ax.arrow(delh_mt*delh+0.3,-1.7+delvar, 0,delar, head_width=hw)
    ax.text(delh_mt*delh+delar+0.3,+lshift-1.7+delvar, r"$x$", fontsize=FontSize)
    ax.text(delh_mt*delh+lshift+0.3,delar-1.7+delvar, r"$z$", fontsize=FontSize)

    ax.arrow(delh_mt*delh,-2.55+delvar, delar,0, head_width=hw)
    ax.arrow(delh_mt*delh,-2.55+delvar, 0,delar, head_width=hw)
    ax.text(delh_mt*delh+delar-lshift-0.2,lshift-2.55+delvar, r"$\phi$", fontsize=FontSize)
    ax.text(delh_mt*delh+lshift,delar-2.55+delvar, r"$z$", fontsize=FontSize)

    ax.arrow(delh_mt*delh,-2.95+delvar, delar,0, head_width=hw)
    ax.arrow(delh_mt*delh,-2.95+delvar, 0,-delar, head_width=hw)
    #ax.text(delh_mt*delh+delar,+lshift-2.95+delvar, r"$\phi$", fontsize=FontSize)
    ax.text(delh_mt*delh+lshift,-delar-2.95+delvar, r"$\rho$", fontsize=FontSize)


    #cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")),ax=ax,ticks=[0,np.pi/6,np.pi/3,np.pi/2])
    ax.set_aspect('equal', 'box')
    ax.set_aspect("equal")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)

    #plt.tight_layout(pad=0)
    plt.subplots_adjust(left=0.05,bottom=0.,right=1,top=1)
    #plt.tight_layout(pad=-1)
    #plt.show()
    plt.savefig("figures/diagram_walls.pdf", format="pdf")



def get_energy_color_map_data():
    pass

def energy_color_map_plot(LineWidth, FontSize, LabelSize):
    pass
    # can use diverging color maps
    # e.g. 'coolwarm'

def get_del_Ftot_data():
    K = 1
    R = 1
    folder = "../data/pydata/Mar20_2023"
    Cs = np.arange(0.2,8.1,0.2)
    optFtot_m0, optalpha_m0, optgamma_m0 = [], [], []
    optFtot_m2, optalpha_m2, optgamma_m2 = [], [], []
    mi,mj = 0, 2
    for C in Cs:
        filename_m0 = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C,mi, R)
        filename_m2 = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C,mj, R)

        qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m0, skiprows=1, delimiter=",", unpack=True)
        optFtot_m0.append(optFtot)
        optalpha_m0.append(optalpha)
        optgamma_m0.append(optgamma)

        qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m2, skiprows=1, delimiter=",", unpack=True)
        optFtot_m2.append(optFtot)
        optalpha_m2.append(optalpha)
        optgamma_m2.append(optgamma)
    return(qs,Cs, optFtot_m0, optalpha_m0, optgamma_m0,optFtot_m2, optalpha_m2, optgamma_m2)

def del_Ftot_phase_Ks_qs_plot(LineWidth, FontSize, LabelSize):
    print("aka: phase diagram for interpolation model")
    qs,Cs,optFtot_m0, optalpha_m0, optgamma_m0,optFtot_m2, optalpha_m2, optgamma_m2 = get_del_Ftot_data()
    mi,mj=0,2
    #qmesh,Kmesh = np.meshgrid(qs,Ks)
    qmesh,Kmesh = np.meshgrid(qs,Cs) # testing CR^2/K
    optFtot_diff = np.array(optFtot_m2) - np.array(optFtot_m0)
    absFdiffmax = np.max(np.abs(optFtot_diff))
    print("absFdiffmax=",absFdiffmax)
    ppi = 72
    plt.figure()
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.65))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axFdiff = plt.subplot2grid((2, 3), (0, 0),rowspan=2, colspan=2)
    axalpha = plt.subplot2grid((2, 3), (0, 2),rowspan=1)
    axgamma = plt.subplot2grid((2, 3), (1, 2),rowspan=1,sharex=axalpha)
    msize = 4

    cmap = cm.get_cmap("bwr")
    #cf = axs[0].pcolormesh(qmesh,Kmesh,optFtot_diff, shading="nearest", cmap = cmap,
                            #norm = colors.DivergingNorm(vmin=np.min(optFtot_diff), vcenter = 0,vmax = np.max(optFtot_diff) ))
    Fdiffmin,Fdiffmax = np.min(optFtot_diff),np.max(optFtot_diff)
    print(Fdiffmin,Fdiffmax)
    #Fdiffmin,Fdiffmax = round(Fdiffmin,0),round(Fdiffmax,2)
    levels = [Fdiffmin,Fdiffmin*3/4,Fdiffmin*2/4,Fdiffmin/4,0,Fdiffmax/4,Fdiffmax*2/4,Fdiffmax*3/4,Fdiffmax]
    print("qmesh.shape(),Kmesh.shape,optFtot.shape()",np.shape(qmesh),np.shape(Kmesh),np.shape(optFtot_diff))
    #cf = axFdiff.pcolormesh(qs,Cs,optFtot_diff,shading="gouraud",cmap=cmap)
    #cf = axFdiff.pcolormesh(qmesh,Kmesh,optFtot_diff,shading="gouraud",cmap=cmap, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))
    #cf = axFdiff.pcolormesh(qmesh,Kmesh,optFtot_diff, shading="nearest",cmap=cmap, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))
    qR = np.linspace(0,1.7,20)
    axFdiff.plot(qR,2*(1+qR**2),linestyle="-",color="gray",label=r"$2(1+q^2R^2)$", zorder=1) # tom's theory
    levels = [Fdiffmin,0,Fdiffmax]
    cfc = axFdiff.contourf(qmesh,Kmesh,optFtot_diff,levels,cmap=cmap,alpha=0.75, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))
    cs = axFdiff.contour(cfc,levels=[0],colors=("k"),linestyles="--",linewidths=(LineWidth))


    axFdiff.legend(loc="center left",ncol=1,columnspacing=0.1,handlelength=0.5,handletextpad=0.5,labelspacing=0.1,markerscale=1,frameon=False,fontsize=LabelSize)
    '''
    divider = make_axes_locatable(axFdiff)
    cax = divider.append_axes("top", size="3%", pad=0.02)
    #cbar=plt.colorbar(cf, cax = cax,ax=axFdiff,orientation="vertical",ticks=levels)
    cbar=plt.colorbar(cf, cax = cax,ax=axFdiff,orientation="horizontal",ticks=levels)
    #cbar.add_lines(cs)
    #cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    cbar.ax.set_xticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,1)),"",str(round(Fdiffmax,1))],fontsize=LabelSize)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(direction="in")
    '''

    #cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    axFdiff.set_xlabel(r"$qR$",fontsize=LabelSize)
    axFdiff.set_ylabel(r"$CR^2/K$",fontsize=LabelSize)
    axFdiff.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=True, labelsize=LabelSize)
    #cbar.ax.set_title(r"$\Delta E'R^2/K$",fontsize=LabelSize)
    axFdiff.xaxis.set_major_locator(MultipleLocator(1))
    axFdiff.xaxis.set_minor_locator(MultipleLocator(0.2))
    axFdiff.yaxis.set_major_locator(MultipleLocator(1))
    axFdiff.yaxis.set_minor_locator(MultipleLocator(0.5))

    if(mi==0 and mj==2):
        axFdiff.text(0.4,6,"Smectic-A\n" + r"$(m=%d)$"%mi,fontsize=LabelSize,multialignment='center')
        axFdiff.text(1.0,2,"Cholesteric\n" + r"$(m=%d)$"%mj,fontsize=LabelSize,multialignment='center')

    x1, y1 = 0.8, 0.05
    axFdiff.text(x1, y1, r"(a)", fontsize=FontSize, transform=axFdiff.transAxes)


    # plot tan alpha
    nlegend = 4
    gap=len(optalpha_m2)//nlegend
    if not gap:
        gap=1
    print("gap",gap)
    ccolors = ["red", "blue", "purple" , "black"]
    c = 0
    for i in range(len(optalpha_m2))[::gap]:
        select = optFtot_diff[i]<0
        #select = optFtot_diff[i]<100000
        axalpha.plot(qs[select][:],np.tan(optalpha_m2[i][select][:]),"-",ms=msize,mfc="None",lw=LineWidth,label="%.1f"%Cs[i],color=ccolors[c])
        #print("K[i]=",Ks[i])
        c=c+1
        print("optalpha_m2[i][select]=",optalpha_m2[i][select])
    #axalpha.set_xlabel(r"$qR$",fontsize=LabelSize)
    axalpha.legend(loc="upper left",ncol=1,columnspacing=0.1,handlelength=0.5,handletextpad=0.5,labelspacing=0.1,markerscale=1,frameon=False,fontsize=LabelSize)
    axalpha.set_ylabel(r"$\tan\alpha$",fontsize=LabelSize)
    axalpha.yaxis.set_label_position("right")
    axalpha.yaxis.tick_right()
    axalpha.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=False, labelleft=False, labelsize=LabelSize)
    axalpha.xaxis.set_major_locator(MultipleLocator(1))
    axalpha.xaxis.set_minor_locator(MultipleLocator(0.5))
    axalpha.yaxis.set_major_locator(MultipleLocator(0.5))
    axalpha.yaxis.set_minor_locator(MultipleLocator(0.25))
    x1, y1 = 0.7, 0.1
    axalpha.text(x1, y1, r"(b)", fontsize=FontSize, transform=axalpha.transAxes)


    # plot gamma
    c = 0
    for i in range(len(optalpha_m2))[::gap]:
        select = optFtot_diff[i]<0
        axgamma.plot(qs[select],optgamma_m2[i][select],"-",ms=msize,mfc="None",lw=LineWidth,color=ccolors[c])
        c+=1
    axgamma.set_xlabel(r"$qR$",fontsize=LabelSize)
    axgamma.set_ylabel(r"$\gamma$",fontsize=LabelSize)
    axgamma.yaxis.set_label_position("right")
    axgamma.yaxis.tick_right()
    axgamma.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=False, labelsize=LabelSize)

    axgamma.xaxis.set_major_locator(MultipleLocator(1))
    axgamma.xaxis.set_minor_locator(MultipleLocator(0.5))

    axgamma.yaxis.set_major_locator(MultipleLocator(0.1))
    axgamma.yaxis.set_minor_locator(MultipleLocator(0.05))

    x1, y1 = 0.7, 0.1
    axgamma.text(x1, y1, r"(c)", fontsize=FontSize, transform=axgamma.transAxes)


    plt.tight_layout(pad=0.01)
    plt.savefig("figures/two_mod_del_E.pdf",format="pdf")
    #plt.show()
    plt.close()

def get_del_Ftot_data_3phase():
    K = 1
    R = 1
    folder = "../data/pydata/Jul20_2023"
    Cs = np.arange(0.0,8.1,0.1)
    optFtot_m0, optalpha_m0, optgamma_m0 = [], [], []
    optFtot_m2, optalpha_m2, optgamma_m2 = [], [], []
    optFtot_m3, optalpha_m3, optgamma_m3 = [], [], []
    mi,mj,mk = 0, 2, 3
    for C in Cs:
        filename_m0 = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C,mi, R)
        filename_m2 = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C,mj, R)
        filename_m3 = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C,mk, R)
        mrow = 201
        jrow = 1
        qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m0, skiprows=1, max_rows=mrow, delimiter=",", unpack=True)
        qs,optFtot,optalpha,optgamma = qs[::jrow],optFtot[::jrow],optalpha[::jrow],optgamma[::jrow]
        optFtot_m0.append(optFtot)
        optalpha_m0.append(optalpha)
        optgamma_m0.append(optgamma)

        qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m2, skiprows=1,  max_rows=mrow, delimiter=",", unpack=True)
        qs,optFtot,optalpha,optgamma = qs[::jrow],optFtot[::jrow],optalpha[::jrow],optgamma[::jrow]
        optFtot_m2.append(optFtot)
        optalpha_m2.append(optalpha)
        optgamma_m2.append(optgamma)

        qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m3, skiprows=1,  max_rows=mrow, delimiter=",", unpack=True)
        print("len(qs)",len(qs))
        qs,optFtot,optalpha,optgamma = qs[::jrow],optFtot[::jrow],optalpha[::jrow],optgamma[::jrow]
        print("len(qs)",len(qs))
        optFtot_m3.append(optFtot)
        optalpha_m3.append(optalpha)
        optgamma_m3.append(optgamma)

    return(qs,Cs, optFtot_m0, optalpha_m0, optgamma_m0, optFtot_m2, optalpha_m2, optgamma_m2, optFtot_m3, optalpha_m3, optgamma_m3)

def m_of_minFtot_3phase_data_get(qs,Cs, optFtot_m0, optFtot_m2, optFtot_m3):
    qmesh,Cmesh = np.meshgrid(qs,Cs) # testing CR^2/K
    allFtot = np.array([optFtot_m0, optFtot_m2, optFtot_m3])
    #Ftot_min = np.min(allFtot,axes)
    allFtot = np.moveaxis(allFtot,[0,1,2],[2,0,1])
    m_of_minFtot = np.argmin(allFtot,axis=2)
    print(m_of_minFtot)
    return m_of_minFtot

def transition_line_of_m_data_get(qs,Cs,m_of_minFtot):
    qm0,Cm0 = [],[]
    qm2,Cm2 = [],[]
    for j in range(len(m_of_minFtot[0])):
        # of different q
        for i in range(len(m_of_minFtot)-1):
            # of diferent C
            if m_of_minFtot[i][j] != 0 and m_of_minFtot[i+1][j] == 0:
                # m0 to other transition
                #qm0.append((qs[j]+qs[j+1])/2)
                #Cm0.append(Cs[i])
                qm0.append(qs[j])
                Cm0.append((Cs[i]+Cs[i+1])/2)

            '''
            elif i<len(m_of_minFtot)-1 and m_of_minFtot[i][j] != 0 and m_of_minFtot[i+1][j] == 0:
                qm0.append(qs[j])
                Cm0.append(Cs[i])
            '''

            if m_of_minFtot[i][j] == 1 and m_of_minFtot[i+1][j] == 2:
                # m2 to m3 transition
                print("m2 to m3: ",qs[j],Cs[i])
                #qm2.append((qs[j]+qs[j+1])/2)
                #Cm2.append(Cs[i])
                qm2.append(qs[j])
                Cm2.append((Cs[i]+Cs[i+1])/2)
            '''
            elif i<len(m_of_minFtot)-1 and m_of_minFtot[i][j] == 1 and m_of_minFtot[i+1][j] == 2:
                print("m2 to m3: ",qs[j],Cs[i])
                qm2.append(qs[j])
                Cm2.append(Cs[i])
            '''
    #qm0 = [qs[0]]+ qm0
    #Cm0 = [Cm0[0]]+ Cm0
    #exlopate the las Cm0
    qm0.append(qm0[-1]+(qm0[-1]-qm0[-2]))
    Cm0.append(Cm0[-1]+(Cm0[-1]-Cm0[-2]))

    print("qm2,Cm2",qm2,Cm2)
    #qm2 = [qs[-1]] + qm2
    #Cm2 = [Cm2[0]]+ Cm2
    print("qm2,Cm2",qm2,Cm2)
    #return (qm0,Cm0,qm2[::-1],Cm2[::-1]) # m2 to 3 transition order is reversed

    # smooth these curves
    Cm0 = savgol_filter(Cm0, 21,3, mode="nearest")
    Cm2 = savgol_filter(Cm2, 21,3, mode="nearest")
    return (qm0,Cm0,qm2,Cm2)




def del_Ftot_phase_Ks_qs_plot_3phase(LineWidth, FontSize, LabelSize):
    print("aka: phase diagram for interpolation model")
    qs,Cs,optFtot_m0, optalpha_m0, optgamma_m0,optFtot_m2, optalpha_m2, optgamma_m2, optFtot_m3, optalpha_m3, optgamma_m3 = get_del_Ftot_data_3phase()
    mi,mj, mk=0,2,3
    qmesh,Cmesh = np.meshgrid(qs,Cs) # testing CR^2/K
    m_of_minFtot = m_of_minFtot_3phase_data_get(qs,Cs, optFtot_m0, optFtot_m2, optFtot_m3) # 0->m=0, 1->m=2, 2->m=3
    qm0,Cm0,qm2,Cm2 = transition_line_of_m_data_get(qs,Cs,m_of_minFtot)
    qm0,Cm0,qm2,Cm2 = np.array(qm0),np.array(Cm0),np.array(qm2),np.array(Cm2)


    ppi = 72
    plt.figure()
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.65))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axFdiff = plt.subplot2grid((2, 3), (0, 0),rowspan=2, colspan=2)
    axalpha = plt.subplot2grid((2, 3), (0, 2),rowspan=1)
    axgamma = plt.subplot2grid((2, 3), (1, 2),rowspan=1,sharex=axalpha)
    msize = 4

    cmap = cm.get_cmap("bwr")

    print(m_of_minFtot)
    # just for confirming the m=0->3 transition
    #axFdiff.pcolormesh(qmesh,Cmesh,m_of_minFtot,shading="auto",cmap=cm.get_cmap("rainbow"))
    #plt.show()


    # plot critical lines
    #axFdiff.plot(qm0,Cm0, "--",color="gray")
    #axFdiff.plot(qm2,Cm2, "--",color="gray")

    alph = 0.6
    # fill m0 section
    axFdiff.fill_between(qm0,Cm0,y2=max(Cs),color="tomato",alpha = alph, edgecolor="None")

    # find tri critial point
    qt,Ct = qm2[0], Cm2[0]
    # fill m2 section
    print("np.concatenate(qm0[qm0<qt],qm2)",np.concatenate((qm0[qm0<qt],qm2)))
    axFdiff.fill_between(np.concatenate((qm0[qm0<qt],qm2)),np.concatenate((Cm0[qm0<qt],Cm2)),y2=min(Cs), color="royalblue",alpha = alph, edgecolor="None")
    #axFdiff.fill_between(qm2,Cm2,y2=min(Cs),color="royalblue",alpha = 0.75)

    # fill m3 section
    axFdiff.fill_betweenx(np.concatenate((Cm2[::-1],Cm0[qm0>=qt])),np.concatenate((qm2[::-1],qm0[qm0>=qt],)),x2=max(qm2), color = "purple",alpha=alph, edgecolor="None")
    #axFdiff.fill_betweenx(Cm2,qm2,x2=max(qm2), color = "purple",alpha=0.75)

    axFdiff.set_xlim(min(qm0),max(qm2))
    axFdiff.set_ylim(0,max(Cs))

    qR = np.linspace(0,1.7,20)
    axFdiff.plot(qR,2*(1+qR**2),linestyle="-",color="k",label=r"$2(1+q^2R^2)$")
    axFdiff.legend(loc="center left",ncol=1,columnspacing=0.1,handlelength=0.5,handletextpad=0.5,labelspacing=0.1,markerscale=1,frameon=False,fontsize=LabelSize)

    '''
    levels = np.array([-0.5,0.5,1.5,2.5])
    cfc = axFdiff.contourf(qmesh,Cmesh,m_of_minFtot,levels, cmap=cmap,alpha=0.75, norm=colors.TwoSlopeNorm(vmin=-1,vcenter=0.5,vmax=3))
    cs = axFdiff.contour(cfc,levels=levels,colors=("k"),linestyles="--",linewidths=(LineWidth))
    qR = np.linspace(0,1.7,20)
    axFdiff.plot(qR,2*(1+qR**2),linestyle="-",color="gray",label=r"$2(1+q^2R^2)$")
    axFdiff.legend(loc="center left",ncol=1,columnspacing=0.1,handlelength=0.5,handletextpad=0.5,labelspacing=0.1,markerscale=1,frameon=False,fontsize=LabelSize)
    '''

    axFdiff.set_xlabel(r"$qR$",fontsize=LabelSize)
    axFdiff.set_ylabel(r"$CR^2/K$",fontsize=LabelSize)
    axFdiff.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=True, labelsize=LabelSize)
    #cbar.ax.set_title(r"$\Delta E'R^2/K$",fontsize=LabelSize)
    axFdiff.xaxis.set_major_locator(MultipleLocator(1))
    axFdiff.xaxis.set_minor_locator(MultipleLocator(0.2))
    axFdiff.yaxis.set_major_locator(MultipleLocator(1))
    axFdiff.yaxis.set_minor_locator(MultipleLocator(0.5))

    if(mi==0 and mj==2 and mk==3):
        #pass
        axFdiff.text(0.4,5,"Smectic-A\n" + r"$(m=%d)$"%mi,fontsize=LabelSize,multialignment='center')
        #axFdiff.text(1.0,2,"Cholesteric\n",fontsize=LabelSize,multialignment='center')
        axFdiff.text(0.4,0.8,"Cholesteric\n" + r"$(m=%d)$"%mj,fontsize=LabelSize,multialignment='center')
        axFdiff.text(1.2,3, "Cholesteric\n" +r"$(m=%d)$"%mk,fontsize=LabelSize,multialignment='center')

    x1, y1 = 0.8, 0.05
    axFdiff.text(x1, y1, r"(a)", fontsize=FontSize, transform=axFdiff.transAxes)

    # to be built for 3 phases case

    # plot tan alpha
    nlegend = 4
    gap=len(optalpha_m2)//nlegend
    if not gap:
        gap=1
    print("gap",gap)
    start = 2
    ccolors = ["gold","red", "blue", "purple", "black"]
    lss = ["--","-"]
    c = 0
    C_to_plot = [1.0,1.5,2.5,4]

    for C in C_to_plot:
        i = np.where(Cs==C)[0][0]
        for m in [2,3]:
            select = m_of_minFtot[i]==m-1
            #select = optFtot_diff[i]<100000
            if m==2:
                optalpha =optalpha_m2[i]
                axalpha.plot(qs[select][:],np.tan(optalpha[select][:]),ls=lss[m-2],ms=msize,mfc="None",lw=LineWidth,label="%.1f"%Cs[i],color=ccolors[c])
            elif m==3:
                optalpha =optalpha_m3[i]
                axalpha.plot(qs[select][:],np.tan(optalpha[select][:]),ls=lss[m-2],ms=msize,mfc="None",lw=LineWidth,color=ccolors[c])

            #print("K[i]=",Ks[i])
        c=c+1
        #print("optalpha_m2[i][select]=",optalpha_m2[i][select])

    #axalpha.set_xlabel(r"$qR$",fontsize=LabelSize)
    axalpha.legend(loc="upper left",ncol=1,columnspacing=0.1,handlelength=0.5,handletextpad=0.25,labelspacing=0.1,markerscale=1,frameon=False,fontsize=LabelSize)
    axalpha.set_ylabel(r"$\tan\alpha$",fontsize=LabelSize)
    axalpha.yaxis.set_label_position("right")
    axalpha.yaxis.tick_right()
    axalpha.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=False, labelleft=False, labelsize=LabelSize)
    axalpha.xaxis.set_major_locator(MultipleLocator(1))
    axalpha.xaxis.set_minor_locator(MultipleLocator(0.5))
    axalpha.yaxis.set_major_locator(MultipleLocator(0.5))
    axalpha.yaxis.set_minor_locator(MultipleLocator(0.25))
    x1, y1 = 0.7, 0.1
    axalpha.text(x1, y1, r"(b)", fontsize=FontSize, transform=axalpha.transAxes)


    # plot gamma
    c = 0
    for C in C_to_plot:
        i = np.where(Cs==C)[0][0]
        for m in [2,3]:
            select = m_of_minFtot[i]==m-1
            if m==2:
                axgamma.plot(qs[select],optgamma_m2[i][select],"--",ms=msize,mfc="None",lw=LineWidth,color=ccolors[c])
            elif m==3:
                axgamma.plot(qs[select],optgamma_m3[i][select],"-",ms=msize,mfc="None",lw=LineWidth,color=ccolors[c])


        c+=1
    axgamma.set_xlabel(r"$qR$",fontsize=LabelSize)
    axgamma.set_ylabel(r"$\gamma$",fontsize=LabelSize)
    axgamma.yaxis.set_label_position("right")
    axgamma.yaxis.tick_right()
    axgamma.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=False, labelsize=LabelSize)

    axgamma.xaxis.set_major_locator(MultipleLocator(1))
    axgamma.xaxis.set_minor_locator(MultipleLocator(0.5))

    axgamma.yaxis.set_major_locator(MultipleLocator(0.1))
    axgamma.yaxis.set_minor_locator(MultipleLocator(0.05))

    x1, y1 = 0.7, 0.1
    axgamma.text(x1, y1, r"(c)", fontsize=FontSize, transform=axgamma.transAxes)


    plt.tight_layout(pad=0.1)
    plt.savefig("figures/two_mod_del_E_3phase.pdf",format="pdf")
    #plt.show()
    plt.close()


