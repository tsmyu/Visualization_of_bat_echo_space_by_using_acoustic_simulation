
import collections
# from minepy import MINE
import math
from scipy.stats import spearmanr
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
import glob
import numpy as np
import matplotlib.colors as mcolor
import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 18
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.labelsize"] = 18
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["legend.loc"] = "best"
plt.rcParams["figure.autolayout"] = True
# import seaborn as sns
# sns.set_style("darkgrid")

dx = 3.0e-4
dt = 5.0e-7
v = 331.5


def read_txt(turn_rate_txt, pulse_gaze_txt, echo_gaze_txt):
    with open(turn_rate_txt) as f:
        turn_rate = [s.strip().split("\t") for s in f.readlines()]
    with open(pulse_gaze_txt) as f:
        pulse_gaze = [s.strip().split("\t") for s in f.readlines()]
    with open(echo_gaze_txt) as f:
        echo_gaze = [s.strip().split("\t") for s in f.readlines()]

    return turn_rate, pulse_gaze, echo_gaze


def calc_corrcoef(turn_rate_pick, gaze_pick):
    turn_rate_arr = np.array(turn_rate_pick)
    gaze_arr = np.array(gaze_pick)
    # turn_rate_arr = turn_rate_arr - np.mean(turn_rate_arr)
    # gaze_arr = gaze_arr - np.mean(gaze_arr)
    score = np.corrcoef(turn_rate_arr, gaze_arr)[0][1]

    return score


def calc_ci(r, n, alpha=0.95):
    z = 0.5*np.log((1+r)/(1-r))
    za = stats.norm.ppf(0.5 + 0.5 * alpha)
    zl = z - za * math.sqrt(1/(n-3))
    zu = z + za * math.sqrt(1/(n-3))
    rhol = (math.exp(2 * zl) - 1) / (math.exp(2 * zl) + 1)
    rhou = (math.exp(2 * zu) - 1) / (math.exp(2 * zu) + 1)
    return rhol, rhou


def zncc(wave1, wave2):
    vec1 = np.array(wave1)
    vec2 = np.array(wave2)
    numer = np.dot(vec1, vec2.T)
    denom = np.sqrt(np.sum(vec1 ** 2)) * np.sqrt(np.sum(vec2 ** 2))
    if denom == 0:
        return 0

    r = numer / denom
    n = len(wave1)
    rhol, rhou = calc_ci(r, n, 0.95)

    return r, rhol, rhou


# def calc_mine(wave1, wave2):
#     mine = MINE()
#     mine.compute_score(wave1, wave2)

#     return mine.mic()


def __make_figure(i, data_name, turn_rate_for_scatter,
                  pulse_gaze_for_scatter, echo_gaze_for_scatter, gaze_tim_for_scatter):

    diff_gaze = np.array(echo_gaze_for_scatter) - \
        np.array(pulse_gaze_for_scatter)
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    colorlist = ["r", "g", "b", "c", "m", "y", "k", "w"]
    for idx in range(len(gaze_tim_for_scatter)):
        ax1.plot([diff_gaze[idx], diff_gaze[idx]], [turn_rate_for_scatter[idx],
                                                    turn_rate_for_scatter[idx]], [0, gaze_tim_for_scatter[idx]], color='gray')

    ax1.set_ylim([-90, 90])
    ax1.set_xlim([-150, 150])

    ax1.scatter(diff_gaze, turn_rate_for_scatter,
                gaze_tim_for_scatter, c=colorlist[i])
    ax1.set_xlabel("Turn rate [degree/t]")
    ax1.set_ylabel("difference from pulse and echo direction [degree]")
    # ax1.set_zlabel("Gaze arrival time [s]")
    plt.grid(True, color='b', linestyle='dotted')
    # plt.savefig(f"{data_name}_diff_scatter.png")
    plt.show()


def __make_fig_tim_axis(data_name, turn_rate_one, pulse_gaze_one, echo_gaze_one):
    fig = plt.figure(tight_layout=True)
    ax1 = fig.add_subplot(
        1, 1, 1, xlabel="time [s]", ylabel="turn rate [degree/t]")
    colorlist = ["grey", "k", "g", "y"]
    ax1.set_ylim([-180, 180])
    turn_rate_max_tim = -100
    for turn_rate in turn_rate_one:
        ax1.scatter(float(turn_rate[0]), float(turn_rate[1]), c="k", s=10)
        turn_rate_max_tim = float(turn_rate[0])
    ax2 = ax1.twinx()
    ax2.set_ylabel("gaze pulse [degree]", rotation=270, labelpad=15)
    ax2.set_ylim([-90, 90])
    for pulse_gaze, echo_gaze in zip(pulse_gaze_one, echo_gaze_one):
        # ax2.scatter(float(gaze[0]), float(gaze[1]), c=colorlist[int(gaze[3])])
        if float(pulse_gaze[0]) <= turn_rate_max_tim:
            # ax2.scatter(float(pulse_gaze[0]), float(pulse_gaze[1]), c="r")
            ax2.scatter(float(pulse_gaze[0]), abs(float(
                echo_gaze[1]) - float(pulse_gaze[1])), c=colorlist[int(echo_gaze[3])])
    plt.grid(True, color='b', linestyle='dotted')
    fig.align_labels()
    plt.savefig(f"{data_name}_diff_abs.png")
    # plt.show()


def __make_fig_pulse_echo(data_name, turn_rate_one, pulse_gaze_one, echo_gaze_one):
    turn_rate_max_tim = -100
    for turn_rate in turn_rate_one:
        turn_rate_max_tim = float(turn_rate[0])
    fig = plt.figure(tight_layout=True)
    ax1 = fig.add_subplot(
        1, 1, 1, xlabel="time [s]", ylabel="gaze [degree]")
    colorlist = ["grey", "k", "g", "y"]
    ax1.set_ylim([-90, 90])
    for pulse_gaze, echo_gaze in zip(pulse_gaze_one, echo_gaze_one):
        if float(pulse_gaze[0]) <= turn_rate_max_tim and float(echo_gaze[0]) <= turn_rate_max_tim:
            ax1.scatter(float(pulse_gaze[0]), float(pulse_gaze[1]), c="r")
            ax1.scatter(float(echo_gaze[0]), float(
                echo_gaze[1]), c=colorlist[int(echo_gaze[3])])
    plt.grid(True, color='b', linestyle='dotted')
    fig.align_labels()
    plt.savefig(f"{data_name}_pulse_echo.png")
    # plt.show()


def calc_corr(turn_rate_txt_list, pulse_gaze_txt_list, echo_gaze_txt_list):
    score_pulse = []
    score_pulse_high = []
    score_pulse_low = []
    score_echo = []
    score_echo_high = []
    score_echo_low = []
    turn_rate_list = []
    pulse_gaze_list = []
    echo_gaze_list = []
    data_name_list = []
    turn_rate_for_scatter = []
    gaze_for_scatter = []
    tau_list = np.linspace(-0.1, 0.6, round(0.7/0.01+1))
    # tau_list = [0]
    for turn_rate_txt, pulse_gaze_txt, echo_gaze_txt in zip(turn_rate_txt_list, pulse_gaze_txt_list, echo_gaze_txt_list):
        data_name = os.path.splitext(os.path.basename(turn_rate_txt))[
            0].split("_")[-1]
        turn_rate, pulse_gaze, echo_gaze = read_txt(
            turn_rate_txt, pulse_gaze_txt, echo_gaze_txt)
        turn_rate_list.append(turn_rate)
        pulse_gaze_list.append(pulse_gaze)
        echo_gaze_list.append(echo_gaze)
        data_name_list.append(data_name)
    # fig = plt.figure()
    # ax1 = fig.add_subplot(1, 1, 1, projection='3d')

    for tau in tau_list:
        turn_rate_pick = []
        pulse_gaze_pick = []
        echo_gaze_pick = []
        turn_rate_with_tim = []
        gaze_with_tim = []
        pulse_tim = []
        echo_tim = []
        for idx, (data_name, turn_rate_one, pulse_gaze_one, echo_gaze_one) in enumerate(zip(data_name_list, turn_rate_list, pulse_gaze_list, echo_gaze_list)):
            turn_rate_for_scatter = []
            pulse_gaze_for_scatter = []
            echo_gaze_for_scatter = []
            gaze_tim_for_scatter = []
            # if tau == -0.5:
            #     # __make_fig_tim_axis(data_name, turn_rate_one,
            #     #                     pulse_gaze_one, echo_gaze_one)
            #     __make_fig_pulse_echo(data_name, turn_rate_one,
            #                           pulse_gaze_one, echo_gaze_one)
            #     plt.close()
            for pulse_gaze, echo_gaze in zip(pulse_gaze_one, echo_gaze_one):
                for turn_rate in turn_rate_one:
                    if float(turn_rate[0])-tau - float(pulse_gaze[0]) <= 0.001 and float(turn_rate[0])-tau - float(pulse_gaze[0]) >= -0.001:
                        turn_rate_pick.append(float(turn_rate[1]))
                        # gaze_pick.append(
                        #     float(echo_gaze[1]) - float(pulse_gaze[1]))
                        pulse_gaze_pick.append(float(pulse_gaze[1]))
                        echo_gaze_pick.append(float(echo_gaze[1]))
                        # for make figure
                        turn_rate_with_tim.append(
                            (float(turn_rate[0]), float(turn_rate[1])))
                        gaze_with_tim.append(
                            (float(pulse_gaze[0]), float(pulse_gaze[1])))
                        if tau == 0:
                            turn_rate_for_scatter.append(float(turn_rate[1]))
                            pulse_gaze_for_scatter.append(float(pulse_gaze[1]))
                            echo_gaze_for_scatter.append(float(echo_gaze[1]))
                            gaze_tim_for_scatter.append(float(echo_gaze[2]))
                            # pulse and echo time
                            pulse_tim.append(float(pulse_gaze[0]))
                            echo_tim.append(float(echo_gaze[0]))
                        break
            if tau == 0:
                print("time difference:", np.mean(np.array(pulse_tim[3:]) - np.array(echo_tim[:-3])))
                # __make_figure(idx, data_name, turn_rate_for_scatter,
                #            pulse_gaze_for_scatter, echo_gaze_for_scatter, gaze_tim_for_scatter)
        r_pulse, rhol_pulse, rhou_pulse = zncc(turn_rate_pick, pulse_gaze_pick)
        r_echo, rhol_echo, rhou_echo = zncc(turn_rate_pick, echo_gaze_pick)


        # r = calc_mine(turn_rate_pick, gaze_pick)
        # turn rate and pulse gaze
        score_pulse.append(r_pulse)
        score_pulse_low.append(rhol_pulse)
        score_pulse_high.append(rhou_pulse)
        # turn rate and echo gaze
        score_echo.append(r_echo)
        score_echo_low.append(rhol_echo)
        score_echo_high.append(rhou_echo)
    # plt.savefig(f"diff_scatter.png")
    plt.close()
    # __make_figure(turn_rate_for_scatter, gaze_for_scatter)
    plt.plot(tau_list, score_pulse, c="r", lw=3)
    plt.plot(tau_list, score_pulse_low, c="r", lw=2, linestyle="dotted")
    plt.plot(tau_list, score_pulse_high, c="r", lw=2, linestyle="dotted")
    plt.plot(tau_list, score_echo, c="b", lw=3)
    plt.plot(tau_list, score_echo_low, c="b", lw=2, linestyle="dotted")
    plt.plot(tau_list, score_echo_high, c="b", lw=2, linestyle="dotted")
    plt.ylim((-1, 1))
    plt.xlim((-0.1, 0.6))
    plt.xlabel(r'$\tau$ [s]')
    plt.ylabel('Correlation coefficient')
    plt.grid(True, color='b', linestyle='dotted')
    plt.show()
    print(f"max pulse index:{tau_list[np.argmax(score_pulse)]}")
    print(f"max echo index:{tau_list[np.argmax(score_echo)]}")


def main():
    argvs = sys.argv
    if len(argvs) < 1:
        print(
            f"Usage: python {argvs[0]} [turn rate folder] [pulse gaze folder] [echo gaze folder]")
        exit()

    turn_rate_list = sorted(glob.glob(f"{argvs[1]}/*1st.*"))
    pulse_gaze_list = sorted(glob.glob(f"{argvs[2]}/*1st.*"))
    echo_gaze_list = sorted(glob.glob(f"{argvs[3]}/*1st.*"))
    calc_corr(turn_rate_list, pulse_gaze_list, echo_gaze_list)


if __name__ == "__main__":
    main()
