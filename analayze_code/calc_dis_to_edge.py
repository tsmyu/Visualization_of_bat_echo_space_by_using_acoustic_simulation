import os
import sys
import glob
from natsort import natsorted
import pandas as pd
import csv
import math
import itertools
import matplotlib.pyplot as plt
# import seaborn as sns
import pandas as pd
import analyzer
import numpy as np
import sympy as sym
from sympy.plotting import plot
from scipy import stats
sym.init_printing(use_unicode=True)

dx = 3.0e-4
dt = 5.0e-7
v = 331.5


def read_csv_data(target_csv):
    emit_data_line = 1000
    echo_data_line_p = 1000
    echo_power_line = 1000
    echo_points_time_r_l_line = 1000
    echo_points_power_r_l_line = 1000
    with open(target_csv, "r", newline='') as f:
        data_map = csv.reader(f)
        for idx, data in enumerate(data_map):
            if idx == emit_data_line:
                if len(data) > 0:
                    emit_data = eval(data[0])
                else:
                    emit_data = ['(-100, -100)']
            if idx == echo_data_line_p:
                if len(data) > 0:
                    echo_data = data
                else:
                    echo_data = ['(-100, -100)']
            if idx == echo_power_line:
                if len(data) > 0:
                    echo_power = data
                else:
                    echo_power = ['0']
            if idx == echo_points_time_r_l_line:
                if len(data) > 0:
                    echo_r_l_tim = data
                else:
                    echo_r_l_tim = []
            if idx == echo_points_power_r_l_line:
                if len(data) > 0:
                    echo_r_l_power = data
                else:
                    echo_r_l_power = []
            if len(data) > 0:
                if data[0] == "emit_points":
                    emit_data_line = idx + 1
                    continue
                elif data[0] == "echo_points":
                    echo_data_line_p = idx + 1
                    continue
                elif data[0] == "echo_points_power":
                    echo_power_line = idx + 1
                    continue
                elif data[0] == "echo_points_time_r_l":
                    echo_points_time_r_l_line = idx + 1
                    continue
                elif data[0] == "echo_points_power_r_l":
                    echo_points_power_r_l_line = idx + 1
                    continue

    return emit_data, echo_data, echo_power, echo_r_l_tim, echo_r_l_power


def read_echo_data(target_csv):
    print("target_csv:", os.path.basename(target_csv))
    emit_data, echo_data, echo_power, echo_r_l_tim, echo_r_l_power = read_csv_data(
        target_csv)
    tim = []
    angle = []
    angle_to_edge_positions = []
    angle_to_edge_times = []
    power = []
    echo_point_list = []

    for (tims, echo_point, powers) in zip(echo_r_l_tim, echo_data, echo_r_l_power):
        tims = eval(tims)
        echo_points = eval(echo_point)
        powers = eval(powers)
        if 0 <= echo_points[0] <= 15000 and 0 <= echo_points[1] <= 5000:
            tim.append((float(tims[0]) + float(tims[1])) / 2)
            # TODO: how to calc degree for echo points
            # angle_tmp = math.degrees(math.atan(
            #     (int(echo_point[0]) - int(emit_data[0])) / (int(echo_point[1] - int(emit_data[1])))))
            angle_tmp = -math.atan((int(echo_points[1]) - int(emit_data[1])) / (
                int(echo_points[0] - int(emit_data[0]))))*180/math.pi
            angle.append(angle_tmp)
            power.append((float(powers[0]) + float(powers[1])) / 2)
            echo_point_list.append(echo_points)
    # calc angle to three edge of obstacle
    if emit_data != ['(-100, -100)']:
        if int(emit_data[0]) <= 5000:
            edge_positions = [((5000, 2800), (5000, 5000)), ((
                8333, 2200), (8333, 0)), ((11666, 2800), (11666, 5000))]
        elif 5000 < int(emit_data[0]) and int(emit_data[0]) <= 8333:
            edge_positions = [((8333, 2200), (8333, 0)),
                              ((11666, 2800), (11666, 5000))]
        elif 8333 < int(emit_data[0]) and int(emit_data[0]) <= 11666:
            edge_positions = [((11666, 2800), (11666, 5000))]
        else:
            edge_positions = []
        if 0 < len(edge_positions):
            for edge in edge_positions:
                # y方向は上下反転しているので-を除く
                if int(edge[0][0]) - int(emit_data[0]) != 0:
                    angle_to_edge_position_1 = math.atan((int(edge[0][1]) - int(5000-emit_data[1])) / (
                        int(edge[0][0]) - int(emit_data[0])))*180/math.pi
                elif int(edge[0][1]) > int(5000-emit_data[1]):
                    angle_to_edge_position_1 = 90
                else:
                    angle_to_edge_position_1 = -90
                if int(edge[1][0]) - int(emit_data[0]) != 0:
                    angle_to_edge_position_2 = math.atan((int(edge[1][1]) - int(5000-emit_data[1])) / (
                        int(edge[1][0]) - int(emit_data[0])))*180/math.pi
                elif int(edge[1][1]) > int(5000-emit_data[1]):
                    angle_to_edge_position_2 = 90
                else:
                    angle_to_edge_position_2 = -90
                angle_to_edge_times_1 = (
                    (int(edge[0][1]) - int(emit_data[1]))**2 + (int(edge[0][0]) - int(emit_data[0]))**2)**0.5*2*dt/v
                angle_to_edge_times_2 = (
                    (int(edge[1][1]) - int(emit_data[1]))**2 + (int(edge[1][0]) - int(emit_data[0]))**2)**0.5*2*dt/v
                angle_to_edge_positions.append(
                    [angle_to_edge_position_1, angle_to_edge_position_2])
                angle_to_edge_times.append(
                    [angle_to_edge_times_1, angle_to_edge_times_2])

    return tim, angle, power, angle_to_edge_positions, angle_to_edge_times, echo_point_list, emit_data


def read_csv_name(csv):
    b_name = os.path.basename(csv)
    folder_name = os.path.splitext(b_name)[0].split("_")[-1]
    print(folder_name)
    if folder_name[-1] == 't':
        flight_num = 0
    elif folder_name[-1] == 'h':
        flight_num = 1
    else:
        print('No flight Num')
        exit
    if folder_name[3] == 'M':
        bat_name = folder_name[3:5]
    else:
        bat_name = folder_name[3]

    return bat_name, flight_num


def calc_bat_number(bat):
    if bat == 'D':
        j = 0 * 7
    elif bat == 'F':
        j = 1 * 7
    elif bat == 'G':
        j = 2 * 7
    elif bat == 'H':
        j = 3 * 7
    elif bat == 'I':
        j = 4 * 7
    elif bat == 'MF':
        j = 5 * 7
    elif bat == 'MH':
        j = 6 * 7
    return j


def read_excel(target_csv, flight_num):
    data = pd.ExcelFile(target_csv)
    input_sheet_name = data.sheet_names
    if flight_num == 0:
        sheet_df = data.parse('first_flight', skiprows=2)
    elif flight_num == 1:
        sheet_df = data.parse('last_flight', skiprows=2)

    return sheet_df


def create_edge_fig(distance_1st, distance_12th):

    data = [sum(distance_1st, []), sum(distance_12th, [])]
    t_value, p_value = stats.ttest_ind(data[0], data[1], equal_var=False)
    print(t_value, p_value)
    print(f'1st num:{len(data[0])}')
    print(f'12th num:{len(data[1])}')
    bottom_1st, up_1st = stats.norm.interval(
        alpha=0.99, loc=np.mean(data[0]), scale=np.std(data[0]))
    print(f'1st: {bottom_1st} < x < {up_1st}')
    bottom_12th, up_12th = stats.norm.interval(
        alpha=0.99, loc=np.mean(data[1]), scale=np.std(data[1]))
    print(f'12th: {bottom_12th} < x < {up_12th}')
    print(f'shapiro 1st:{stats.shapiro(data[0])}')
    print(f'shapiro 12th:{stats.shapiro(data[1])}')
    print(
        f'U test:{stats.mannwhitneyu(data[0], data[1], use_continuity=True, alternative=None)}')
    value_num_1st = []
    value_num_12th = []
    x = []
    for i in range(8):
        value_num_1st.append(np.count_nonzero(
            (i*10 <= np.array(data[0])) & (np.array(data[0]) < (i+1)*10)))
        value_num_12th.append(np.count_nonzero(
            (i*10 <= np.array(data[1])) & (np.array(data[1]) < (i+1)*10)))
        x.append(i)
    print(x)
    print(f'1st dist:{value_num_1st[0]/len(data[0])}')
    print(f'12th dist:{value_num_12th[0]/len(data[1])}')
    # histgram
    x = np.array(x)
    width = 0.2
    # p_1st = plt.bar(x, value_num_1st, width=width, color='lightgreen')
    # p_12th = plt.bar(x+width, value_num_12th, width=width, color='forestgreen')
    # # plt.xticks(x+width/2, ['0', '10', '20',
    # #                        '30', '0.4', '0.5', '0.6', '0.7'])
    # # plt.xticks(x[:3]+width/2, ['0', '10', '20', '30'])
    # plt.xlabel('Distance to inner edges [cm]')
    # plt.ylabel('Echo incidence point number')
    # plt.legend((p_1st[0], p_12th[0]), ('First flight', 'Last flight'))
    # plt.show()
    # normalize
    p_1st = plt.bar(x, np.array(value_num_1st) /
                    len(data[0]), width=width, color='lightgreen')
    p_12th = plt.bar(x+width, np.array(value_num_12th) /
                     len(data[1]), width=width, color='forestgreen')
    plt.xticks(x+width/2, ['0', '10', '20',
                           '30', '40', '50', '60', '70'])
    plt.xlabel('Distance to inner edge [cm]')
    plt.ylabel('Probability')
    plt.legend((p_1st[0], p_12th[0]), ('First flight', 'Last flight'))
    plt.xlim([-0.1, 3.5])
    plt.ylim([0, 1.0])
    plt.show()
    # heights = [np.mean(data[0]), np.mean(data[1])]
    # bars = np.arange(len(heights))
    # std = [np.std(data[0]), np.std(data[1])]
    # fig = plt.figure()
    # ax1 = fig.add_subplot(1, 1, 1)
    # # ax1.violinplot(data, showmedians=True)
    # ax1.boxplot(data, patch_artist=True,
    #             boxprops=dict(facecolor='g',  # boxの塗りつぶし色の設定
    #                           color='black', linewidth=3),  # boxの枠線の設定
    #             medianprops=dict(color='black', linewidth=3),  # 中央値の線の設定
    #             whiskerprops=dict(color='black', linewidth=3),  # ヒゲの線の設定
    #             capprops=dict(color='black', linewidth=5),  # ヒゲの先端の線の設定
    #             labels=['First flight', 'Last flight'], showfliers=False)
    # ax1.set_ylabel('Distance to inner edges [m]')
    # ax1.set_ylim([0, 1.0])
    # plt.show()
    # plt.close(fig)

def dis_to_arr(dis_data_list):
    data_set = []
    for dis_data in dis_data_list:
        bat_name = dis_data[0]
        flight_num = dis_data[1]
        for pulse_dis in dis_data[2]:
            data_set.append([bat_name, pulse_dis[0], pulse_dis[1], flight_num])
    
    return np.array(np.vstack(data_set))


def main(csv_list, all_data_excel):
    distance_edge_1st = []
    distance_edge_12th = []
    flight_dir_1st = []
    pulse_dir_1st = []
    flight_dir_12th = []
    pulse_dir_12th = []
    corr_1st = []
    corr_12th = []
    dis_only_1st = []
    dis_only_12th = []
    for csv in csv_list:
        echo_list = []
        emit_list = []
        bat_name, flight_num = read_csv_name(csv)
        j = calc_bat_number(bat_name)
        sheet_df = read_excel(all_data_excel, flight_num)
        print(f'target bat:{bat_name}-{flight_num}')
        echo_point_list = natsorted(glob.glob("{}/*.csv".format(csv)))
        for idx, echo_data in enumerate(echo_point_list):
            tim, angle, power, angle_to_edge_positions, angle_to_edge_times, echo_point_list, emit_data = read_echo_data(
                echo_data)
            echo_list.append(
                (idx, tim, angle, power, angle_to_edge_positions, angle_to_edge_times, echo_point_list))
            emit_list.append(emit_data)
        bat_one_flight = sheet_df.iloc[:, [
            j+1, j+2, j+3, j+4, j+5, j+6]].values.tolist()
        ana = analyzer.Analyzer(
            emit_list, echo_list, bat_one_flight, bat_name, flight_num)
        distance_edge, dis_only = ana.calc_base_flight_direction()
        if flight_num == 0:
            distance_edge_1st.append((bat_name, flight_num, distance_edge))
            dis_only_1st.append(dis_only)
        else:
            distance_edge_12th.append((bat_name, flight_num, distance_edge))
            dis_only_12th.append(dis_only)

    dis_1st_arr = dis_to_arr(distance_edge_1st)
    dis_12th_arr = dis_to_arr(distance_edge_12th)
    df_1st = pd.DataFrame(dis_1st_arr)
    df_12th = pd.DataFrame(dis_12th_arr)
    df_1st.to_csv('1st.csv')
    df_12th.to_csv('12th.csv')
    create_edge_fig(dis_only_1st, dis_only_12th)


if __name__ == '__main__':
    argvs = sys.argv
    if len(argvs) < 2:
        print(
            f'Usage: python {argvs[0]} [folder of echo csv list] [data excel]')
        exit()

    csv_folder_list = natsorted(glob.glob(f'{argvs[1]}/echo_point_*'))
    all_data = argvs[2]
    main(csv_folder_list, all_data)
