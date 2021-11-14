#!/opt/homebrew/bin/python3
from plot import *
from catelinder_cal import *
from catelinder_plot import *
from unduloid_cal import *
from unduloid_plot import *
from catenoid_cal import *

def main():
    print("python code for numerical calculation of pulling experiment o/n unduloid")
    #A_z1_e_theta1_demo()
    #res = opt_unduloid_E(5.0, 40, 0.1, 35, 300)

    #unduloid_A_l_plot()
    #unduloid_F_plot()
    #unduloid_config_plot(res.x[0], res.x[1], res.x[2])

    #res = opt_catelinder_E(2.0, 10, 10, 20, 300)
    #print(res)
    #catelinder_config_plot(res.x[0],res.x[1],res.x[2])
    #print("catelinder_A",catelinder_A(res.x[0], res.x[1], res.x[2]))
    #catelinder_F_plot()
    catenoid_A_L_plot()

if __name__ == "__main__":
    main()

