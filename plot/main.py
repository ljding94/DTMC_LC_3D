#!/opt/homebrew/bin/python3
from simple_pulling import *
from nematic_wall import *
from config_plot import *
from topo_chiral import *

def main():
    print("(施工中👷‍♀️) plotting publication-quality figures")
    LineWidth, FontSize, LabelSize = 1,9,8
    #force_pull_plot(LineWidth, FontSize, LabelSize)
    #tilt_Kd_plot(LineWidth, FontSize, LabelSize)
    #twist_q_plot(LineWidth, FontSize, LabelSize)
    topo_change_plot(LineWidth, FontSize, LabelSize)

    #catenoid_trinoid_demo_plot(LineWidth, FontSize, LabelSize)

    # for slides
    #catenoid_demo_3d_plot()
if __name__ == "__main__":
    main()
