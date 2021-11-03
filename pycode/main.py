#!/opt/homebrew/bin/python3
from plot import *
from catelinder_cal import *
from catelinder_plot import *

def main():
    print("python code for numerical calculation of pulling experiment on unduloid")
    #A_z1_e_theta1_demo()
    res = opt_bh_E(5, 40, 15, 600)
    print(res)
    print("res.fun",res.fun)
    print("res.x",res.x)

    catelinder_F_plot()
    catelinder_config_plot(res.x[0],res.x[1],res.x[2])
if __name__ == "__main__":
    main()

