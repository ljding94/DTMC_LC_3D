#!/opt/homebrew/bin/python3
from simple_pulling import *
from nematic_wall import *
from twisting_wall import *
from config_plot import *
from topo_chiral import *

def main():

    print("(施工中👷‍♀️) plotting publication-quality figures")
    LineWidth, FontSize, LabelSize = 1,9,8
    #catenoid_demo(LineWidth, FontSize, LabelSize)
    #init_config_demo(LineWidth, FontSize, LabelSize)
    #force_pull_plot(LineWidth, FontSize, LabelSize)
    #tilt_Kd_plot(LineWidth, FontSize, LabelSize)
    #twist_q_plot(LineWidth, FontSize, LabelSize)
    topo_change_plot(LineWidth, FontSize, LabelSize)

if __name__ == "__main__":
    main()


# data in use:
'''
catenoid_demo: "../data/Ne2/Feb_2022/Feb28_2022/State_N300_imod3_Ne2_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_q1.4_Cn4.0.csv"
init_config_demo: init_config folder
force_pull_plot: no LC  and with LC: "../data/Ne2/May13_2022"  15k MC step lam 6
tilt_Kd_plot: "../data/Ne2/May12_2022"  15k MC step lam 6
twist_q_plot: "../data/Ne2/May12_2022"  15k MC step lam 6
topo_change_plot: "../data/Ne2/Feb_2022/Feb28_2022" 15k MC step kar30 lam6 Kd4
'''
