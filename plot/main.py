#!/usr/local/bin/python3
from simple_pulling import *
from nematic_wall import *
from cholesteric_wall import *
from twisting_wall import *
from config_plot import *
from topo_chiral import *
from twomod_plot import *


def main():

    print("(施工中👷‍♀️) plotting publication-quality figures")
    LineWidth, FontSize, LabelSize = 1,9,8
    #config_demo(LineWidth, FontSize, LabelSize)

    #init_config_demo(LineWidth, FontSize, LabelSize)
    #topo_change_plot(LineWidth, FontSize, LabelSize)

    #topo_Hdis_plot(LineWidth, FontSize, LabelSize)

    #cylinder_config_demo(LineWidth, FontSize, LabelSize)
    #force_pull_plot(LineWidth, FontSize, LabelSize)
    #nematic_Kd_plot(LineWidth, FontSize, LabelSize)
    #cholesteric_q_plot(LineWidth, FontSize, LabelSize)

    walls_Cn_lf_vs_Kd_q(LineWidth, FontSize, LabelSize) #smectic_to_walls


    #walls_3phase_diagram(LineWidth, FontSize, LabelSize) #phase_diagram_wall

    #twist_q_plot(LineWidth, FontSize, LabelSize) # 2 to 3 wall transition

    #wall_pitch_q_plot(LineWidth, FontSize, LabelSize) # break in to the following two plot

    #wall_pitch_q_analysis_plot(LineWidth, FontSize, LabelSize)

    #wall_pitch_q_result_plot(LineWidth, FontSize, LabelSize)

    #del_Ftot_phase_Ks_qs_plot(LineWidth, FontSize, LabelSize)

    #del_Ftot_phase_Ks_qs_plot_3phase(LineWidth, FontSize, LabelSize)

    #demo_config_2mod_u_plot(LineWidth, FontSize, LabelSize)

    #two_mod_diagram_plot(LineWidth, FontSize, LabelSize)

    #walls_membrane_shape_slice(LineWidth, FontSize, LabelSize)

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
