from unduloid import *
import matplotlib.pyplot as plt
from scipy import interpolate


def A_z1_e_theta1_demo():
    e,theta1=np.meshgrid(np.linspace(0.1,0.9,100),np.linspace(0.1*np.pi,0.9*np.pi,100))
    A = []
    for i in range(len(e)):
        A.append([])
        for j in range(len(e[0])):
            A[i].append(unduloid_A_unit_a(e[i][j], theta1[i][j]))
    A = np.array(A)
    lf = unduloid_lf_unit_a(e, theta1)
    print("np.shape(A),np.shape(lf)",np.shape(A),np.shape(lf))
    Af=A.flatten()
    lff = lf.flatten()

    fig = plt.figure()
    ax1 = fig.add_subplot(131,projection='3d')
    ax2 = fig.add_subplot(132,projection='3d')
    ax3 = fig.add_subplot(133,projection='3d')

    ax1.plot_surface(e, lf, A)
    ax1.set_zlabel("A")
    ax2.plot_surface(e, lf, theta1)
    ax2.set_zlabel("theta1")

    f = interpolate.interp2d(e,lf,theta1,kind="cubic")

    plt.show()
    plt.close()


