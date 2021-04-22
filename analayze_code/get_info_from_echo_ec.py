import os
import sys
import glob
import numpy as np
from scipy import signal
from scipy.signal import argrelmax
import itertools
import math
import csv
import sympy as sp
from sympy import sqrt, symbols
import re
import matplotlib.pyplot as plt

# simulation parameter
dx = 3.0e-4
dt = 5.0e-7
# distance between nose and ears [m]
distance_ears = 0.02
velocity_air = 331.5
limit_time = distance_ears / velocity_air

# not delete area
# ignore area : 5cm(5)
ignore_time = 100e-3*2/velocity_air
ignore_points = round(ignore_time/dt/10)

# wave for convolution


class DataField:
    def __init__(self, input_data, emit_csv_list):
        # threash
        self.threash = 0.02
        self.base_name = input_data
        file_name = os.path.basename(input_data)
        print(file_name)
        # self.pixel_area = (
        #     int(file_name.split("_")[-3]), int(file_name.split("_")[-2]))
        self.emit_angle = int(
            re.split("[g]", file_name.split("_")[-1].split(".")[0])[-1])
        self.emit_point = (int(os.path.basename(input_data).split("_")[-3]), int(
            os.path.basename(input_data).split("_")[-2]))
        for emit_csv in emit_csv_list:
            emit_name = os.path.basename(emit_csv)
            _emit_angle = int(
                re.split("[g]", emit_name.split("_")[-1].split(".")[0])[-1])
            if self.emit_angle == _emit_angle:
                emit_data_list = np.genfromtxt(emit_csv, usecols=(
                    0, 1, 2, 3), skip_header=1, skip_footer=1, delimiter=",")
                break
        echo_data_list = np.genfromtxt(input_data, usecols=(
            0, 1, 2, 3), skip_header=1, skip_footer=1, delimiter=",")

        self.time_line = echo_data_list[:, 0]
        self.emit_wave = echo_data_list[:, 1]
        right_wave = echo_data_list[:, 3]
        left_wave = echo_data_list[:, 2][:len(right_wave)]
        emit_zero = np.zeros(len(right_wave))
        emit_right = np.concatenate([emit_data_list[:, 3], emit_zero])[
            :len(right_wave)]
        emit_left = np.concatenate([emit_data_list[:, 3], emit_zero])[
            :len(right_wave)]
        right_wave = right_wave - emit_right
        left_wave = echo_data_list[:, 2] - emit_left
        wave = self.__create_waves()
        self.emit_wave = np.convolve(self.emit_wave, wave, mode='full')
        right_wave = np.convolve(right_wave, wave, mode='full')
        left_wave = np.convolve(left_wave, wave, mode='full')
        self.echo_without_emit_wave = [right_wave, left_wave]

        # 取得・計算する内部変数
        self.data_len = 0
        self.right_corr_raw = None
        self.left_corr_raw = None
        self.right_corr = None
        self.left_corr = None
        self.right_echo_time = None
        self.left_echo_time = None
        self.right_echo_power = None
        self.left_echo_power = None
        self.distance_list = []
        self.angle_list = []
        self.power_list = []

    def __create_waves(self):
        f1 = 68e3         # start frequency
        f2 = 50e3         # end frequency
        fs = int(1/dt)        # sampling freq
        d = 2e-3  # duration d*fsで測定したデータの個数
        time = np.linspace(0, 1, fs, endpoint=False)  # array を生成　0~1までfs*10個の点
        time_sig = time[:int(d*fs)]
        # swav = np.sin((2*np.pi*f1*(f1/f2)**(1/d)**time_sig)/np.log((f1/f2)**(1/d)))#exponential型
        FM2 = np.sin(2*np.pi*(f1*time_sig+(f2-f1)*time_sig**2/2/d))  # linear型
        FM1 = np.sin(2*np.pi*(f1*time_sig+(f2/2-f1/2)
                              * time_sig**2/2/d))/100  # linear型
        FM3 = np.sin(2*np.pi*(f1*time_sig+(3*f2/2-3*f1/2)
                              * time_sig**2/2/d))/100  # linear型
        CF2 = np.sin(2*np.pi*(f1*time_sig))
        CF1 = np.sin(2*np.pi*(f1/2*time_sig))/100
        CF3 = np.sin(2*np.pi*(3*f1/2*time_sig))/100
        # FM or CF
        CF = CF1 + CF2 + CF3
        FM = FM1 + FM2 + FM3
        # waves = np.concatenate([CF, FM])
        waves = FM
        # waves = np.fliplr([waves])[0]

        return waves

    def preprocessing(self):
        # correlation
        corr_wave = self.__get_correlation()
        # envelope
        corr_wave_envelop = self.__get_envelope(corr_wave)
        # notmalize
        corr_wave_envelop_norm = self.__normalization(
            corr_wave_envelop, "emit")

        # set length of data
        self.data_len = min(len(self.echo_without_emit_wave[0]), len(
            corr_wave[1]), len(self.time_line))

        self.right_corr_raw = corr_wave_envelop[1][:self.data_len]
        self.left_corr_raw = corr_wave_envelop[2][:self.data_len]

        self.right_corr = corr_wave_envelop_norm[1][:self.data_len]
        self.left_corr = corr_wave_envelop_norm[2][:self.data_len]
        # import matplotlib as mpl
        # mpl.use('tkagg')
        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax_right = fig.add_subplot(211)
        # ax_left = fig.add_subplot(212)
        # ax_right.plot(corr_wave[1])
        # ax_left.plot(corr_wave[2])
        # ax_right.plot(corr_wave_envelop[1])
        # ax_left.plot(corr_wave_envelop[2])
        # plt.savefig(f"{os.path.basename(self.base_name)}.png")

    def __get_correlation(self):
        """
        calc correlation value
        """
        corr_emit = np.correlate(self.emit_wave, self.emit_wave, "full")[
            len(self.emit_wave):]
        corr_right = np.correlate(self.echo_without_emit_wave[0], self.emit_wave, "full")[
            len(self.emit_wave):]
        corr_left = np.correlate(self.echo_without_emit_wave[1], self.emit_wave, "full")[
            len(self.emit_wave):]

        return [corr_emit, corr_right, corr_left]

    def __get_envelope(self, corr_wave):
        envelope_emit = abs(signal.hilbert(corr_wave[0]))
        envelope_right = abs(signal.hilbert(corr_wave[1]))
        envelope_left = abs(signal.hilbert(corr_wave[2]))

        return [envelope_emit, envelope_right, envelope_left]

    def __normalization(self, wave_list, base="emit"):
        """
        normalize data
        """
        norm_list = []
        if base == "emit":
            for idx, wave in enumerate(wave_list):
                if idx == 0:
                    self.wave_max = max(np.array(wave))/100
                wave = np.array(wave)
                wave = wave - wave.mean()
                wave_n = wave / self.wave_max
                norm_list.append(wave_n)
        elif base == "echo":
            for wave in wave_list:
                wave = np.array(wave)
                wave = wave - wave.mean()
                wave_max = max(wave)
                wave_n = wave / wave_max
                norm_list.append(wave_n)

        return norm_list

    def get_info(self):
        right_echo_point, left_echo_point = self.__get_echo_point()
        self.right_echo_time, self.left_echo_time = self.__get_echo_time(
            right_echo_point, left_echo_point)
        self.right_echo_power, self.left_echo_power = self.__get_echo_power(
            right_echo_point, left_echo_point)
        position_list, r_l_tim_list, power_list, r_l_power_list = self.__get_position()
        self.power_list = power_list
        self.r_l_tim_list = r_l_tim_list
        self.r_l_power_list = r_l_power_list
        self.echo_points = self.__get_echo_points(position_list)
        # tmp
        # # 2020/10/08
        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax1 = fig.add_subplot(211)
        # ax2 = fig.add_subplot(212)
        # ax1.plot(self.right_corr)
        # ax1.scatter(right_echo_point, self.right_echo_power /
        #             self.wave_max, c="r")
        # ax2.plot(self.left_corr)
        # ax2.scatter(left_echo_point, self.left_echo_power/self.wave_max, c="r")
        # path = os.getcwd()
        # s_folder = os.path.basename(self.base_name).split("_")[2]
        # ss_folder = os.path.basename(self.base_name).split("_")[3]
        # if not os.path.exists(f"{path}/{s_folder}/{ss_folder}"):
        #     os.makedirs(f"{path}/{s_folder}/{ss_folder}")
        # file_name = os.path.basename(self.base_name).split(".")[0]
        # plt.savefig(f"{path}/{s_folder}/{ss_folder}/{file_name }.png")

    def __get_echo_point(self):
        # right_echo_point_raw = np.where(self.right_corr >= self.threash)[0]
        # left_echo_point_raw = np.where(self.left_corr >= self.threash)[0]
        right_echo_point_raw = self.right_corr
        left_echo_point_raw = self.left_corr

        # pre_r_point = 0
        # pre_point = 0
        # for r_point in right_echo_point_raw:
        #     if r_point - pre_r_point != 1:
        #         if pre_point == 0:
        #             pre_point = r_point
        #         elif pre_point != 0:
        #             echo_point = round((pre_point+r_point)/2)
        #             right_echo_list.append(int(echo_point))
        #             pre_point = 0
        #     pre_r_point = r_point
        # pre_l_point = 0
        # pre_point = 0
        # for l_point in left_echo_point_raw:
        #     if l_point - pre_l_point != 1:
        #         if pre_point == 0:
        #             pre_point = l_point
        #         elif pre_point != 0:
        #             echo_point = round((pre_point+l_point)/2)
        #             left_echo_list.append(int(echo_point))
        #             pre_point = 0
        #     pre_l_point = l_point
        peaks_right = argrelmax(right_echo_point_raw, order=100)
        peak_power_right = right_echo_point_raw[peaks_right]
        peak_right = peaks_right[0][np.where(
            peak_power_right > self.threash)[0]]
        peaks_left = argrelmax(left_echo_point_raw, order=100)
        peak_power_left = left_echo_point_raw[peaks_left]
        peak_left = peaks_left[0][np.where(peak_power_left > self.threash)[0]]
        # import matplotlib as mpl
        # mpl.use('tkagg')
        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax_right = fig.add_subplot(211)
        # ax_left = fig.add_subplot(212)
        # ax_right.plot(right_echo_point_raw)
        # ax_right.scatter(peak_right, right_echo_point_raw[peak_right], c="r")
        # ax_left.plot(left_echo_point_raw)
        # ax_left.scatter(peak_left, left_echo_point_raw[peak_left], c ="r")
        # plt.savefig(f"{os.path.basename(self.base_name)}.png")

        return peak_right, peak_left

    def __get_echo_time(self, right_points, left_points):

        return self.time_line[right_points], self.time_line[left_points]

    def __get_echo_power(self, right_points, left_points):

        return self.right_corr_raw[right_points], self.left_corr_raw[left_points]

    def __get_position(self):
        if min(len(self.right_echo_time), len(self.left_echo_time)) == 0:
            position_list = [[(-100, -100)]]
            r_l_tim_list = [[(-100, -100)]]
            power_list = []
            r_l_power_list = []
        else:
            position_list, r_l_tim_list, power_list, r_l_power_list = self.__calc_position()

        return position_list, r_l_tim_list, power_list, r_l_power_list

    def __calc_position(self):
        position_list = []
        r_l_tim_list = []
        power_list = []
        r_l_power_list = []
        x = symbols("x", positive=True)
        y = symbols("y", positive=True)
        a1 = symbols("a1", positive=True)
        a2 = symbols("a2", positive=True)
        sp.var('x, y, a1, a2')
        eq1 = (x-(distance_ears / 4)) ** 2 / a1 ** 2 + \
            y ** 2 / (a1 ** 2 - (distance_ears / 4) ** 2) - 1
        eq2 = (x+(distance_ears / 4)) ** 2 / a2 ** 2 + \
            y ** 2 / (a2 ** 2 - (distance_ears / 4) ** 2) - 1
        ans = sp.solve([eq1, eq2], [x, y])
        for right_tim, right_power in zip(self.right_echo_time, self.right_echo_power):
            tmp_position_list = []
            tmp_r_l_tim_list = []
            tmp_power_list = []
            tmp_r_l_power_list = []
            for left_tim, left_power in zip(self.left_echo_time, self.left_echo_power):
                if abs(right_tim - left_tim) < limit_time:
                    a_1 = right_tim * velocity_air / 2
                    a_2 = left_tim * velocity_air / 2
                    for i in range(len(ans)):
                        x0 = ans[i][0].subs([(a1, a_1), (a2, a_2)])
                        y0 = ans[i][1].subs([(a1, a_1), (a2, a_2)])
                        if x0.is_real and y0.is_real:
                            if y0 >= 0:
                                x = round(x0 / dx)
                                y = round(y0 / dx)
                    # x = round(ans[1][0].subs([(a1, a_1), (a2, a_2)]) / dx)
                    # y = round(ans[1][1].subs([(a1, a_1), (a2, a_2)]) / dx)
                    tmp_position_list.append([x, y])
                    tmp_r_l_tim_list.append([right_tim, left_tim])
                    tmp_power_list.append(((right_power + left_power)))
                    tmp_r_l_power_list.append([right_power, left_power])
            position_list.append(tmp_position_list)
            r_l_tim_list.append(tmp_r_l_tim_list)
            power_list.append(tmp_power_list)
            r_l_power_list.append(tmp_r_l_power_list)
        r_l_tim_l = list(itertools.chain(*r_l_tim_list))
        r_l_tim_list = sorted(r_l_tim_l, key=r_l_tim_l.index)
        power_l = list(itertools.chain(*power_list))
        power_list = sorted(power_l, key=power_l.index)
        r_l_power_l = list(itertools.chain(*r_l_power_list))
        r_l_power_list = sorted(r_l_power_l, key=r_l_power_l.index)
        return position_list, r_l_tim_list, power_list, r_l_power_list

    def __get_echo_points(self, position_list):
        pixel_points = []
        theta = math.radians(-self.emit_angle)
        rotation_matrix = np.array(
            [np.cos(theta), -np.sin(theta), np.sin(theta), np.cos(theta)]).reshape(2, 2)
        print(position_list)
        for position_onepulse in position_list:
            for position in position_onepulse:
                print(position)
                if position == []:
                    pass
                else:
                    pos_arr = np.array(position).reshape(2, 1)
                    point_arr = np.dot(rotation_matrix, pos_arr)
                    pixel_points.append((round(
                        point_arr[0][0] + self.emit_point[0]), round(self.emit_point[1] - point_arr[1][0])))

        return pixel_points

    def STFT(self, wav_data):
        """
        STFT関数
        input: wav_data
        output: spectrogram
        """
        fs = int(1/dt)
        f, t, spectrogram = signal.stft(
            wav_data.reshape(-1), fs=fs, nperseg=2048, noverlap=2000)

        return f, t, spectrogram

    def save(self):
        print(self.power_list)
        with open("./echo_point_{}/echo_point_{}.csv".format(os.path.basename(os.path.dirname(self.base_name)),
                                                             os.path.splitext(self.base_name)[0].split("\\")[-1]), "a") as f:
            writer = csv.writer(f, lineterminator='\n')

            writer.writerow(["right_echo_point"])
            writer.writerow(self.right_echo_time)
            writer.writerow(["left_echo_point"])
            writer.writerow(self.left_echo_time)
            writer.writerow(["distance"])
            writer.writerow(self.distance_list)
            writer.writerow(["angle"])
            writer.writerow(self.angle_list)
            writer.writerow(["emit_points"])
            writer.writerow([self.emit_point])
            writer.writerow(["echo_points"])
            writer.writerow(self.echo_points)
            writer.writerow(["echo_points_power"])
            writer.writerow(self.power_list)
            writer.writerow(["echo_points_time_r_l"])
            writer.writerow(self.r_l_tim_list)
            writer.writerow(["echo_points_power_r_l"])
            writer.writerow(self.r_l_power_list)

        tim = self.time_line[:self.data_len][:, np.newaxis]
        echo_right = self.echo_without_emit_wave[0][:self.data_len][:, np.newaxis]
        echo_left = self.echo_without_emit_wave[1][:self.data_len][:, np.newaxis]
        corr_right = self.right_corr[:, np.newaxis]
        corr_left = self.left_corr[:, np.newaxis]
        data_for_csv = np.c_[tim, echo_right, echo_left, corr_right, corr_left]
        # for figure
        plt.rcParams["font.size"] = 32
        plt.rcParams["axes.labelsize"] = 32
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax_right = fig.add_subplot(
            111, xlabel="Time [ms]", ylabel="Sound pressure", xticks=[0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0], yticks=[-0.01, 0, 0.01])
        ax_right.plot(self.time_line[:self.data_len]*1000, echo_right)
        plt.savefig(f"wave/right_wave_{os.path.basename(self.base_name)}.png")
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax_stft_right = fig.add_subplot(
            111, xlabel="Time [ms]", ylabel="Frequency [kHz]")
        f_axis, t_axis, spectro = self.STFT(echo_right)
        print(f'f_axis shape:{f_axis.shape}')
        ax_stft_right.pcolormesh(t_axis, f_axis, np.abs(spectro))
        plt.ylim([0, 150e3])
        ax_stft_right.set_yticklabels([0, 25, 50, 75, 100, 125, 150])
        ax_stft_right.set_xticklabels(
            [0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0])
        plt.savefig(
            f"wave/stft_right_wave_{os.path.basename(self.base_name)}.png")
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        # ax_left = fig.add_subplot(
        #     111, xlabel="time [s]", ylabel="sound pressure", yticks=[-0.01, 0, 0.01])
        ax_left = fig.add_subplot(
            111, xlabel="Time [ms]", ylabel="Sound pressure", xticks=[0, 5.0, 10.0, 15.0, 20.0], yticks=[-0.01, 0, 0.01])
        ax_left.plot(self.time_line[:self.data_len]*1000, echo_left)
        plt.savefig(f"wave/left_wave_{os.path.basename(self.base_name)}.png")
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax_stft_left = fig.add_subplot(
            111, xlabel="Time [ms]", ylabel="Frequency [kHz]")
        f_axis, t_axis, spectro = self.STFT(echo_left)
        print(f'f_axis shape:{f_axis.shape}')
        ax_stft_left.pcolormesh(t_axis, f_axis, np.abs(spectro))
        plt.ylim([0, 150e3])
        ax_stft_left.set_yticklabels([0, 25, 50, 75, 100, 125, 150])
        ax_stft_left.set_xticklabels(
            [0, 5.0, 10.0, 15.0, 20.0])
        plt.savefig(
            f"wave/stft_left_wave_{os.path.basename(self.base_name)}.png")
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax_emit = fig.add_subplot(
            111, xlabel="Time [ms]", ylabel="Sound pressure", xticks=[0, 5.0, 10.0, 15.0, 20.0])
        ax_emit.plot(self.time_line*1000,
                     self.emit_wave[:len(self.time_line)])
        plt.savefig(f"emit/emit_{os.path.basename(self.base_name)}.png")
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax_stft_emit = fig.add_subplot(
            111, xlabel="Time [s]", ylabel="Frequency [kHz]")
        f_axis, t_axis, spectro = self.STFT(
            self.emit_wave[:len(self.time_line)])
        print(f'f_axis shape:{f_axis.shape}')
        ax_stft_emit.pcolormesh(t_axis, f_axis, np.abs(spectro))
        plt.ylim([0, 150e3])
        ax_stft_emit.set_yticklabels([0, 25, 50, 75, 100, 125, 150])
        ax_stft_emit.set_xticks(
            [0, 5e-3, 10e-3, 15e-3, 20e-3])
        ax_stft_emit.set_xticklabels(
            [0, 5.0, 10.0, 15.0, 20.0])
        plt.savefig(f"emit/stft_emit_{os.path.basename(self.base_name)}.png")
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax_right_tau = fig.add_subplot(
            111, xlabel=r'$\tau$ [ms]', ylabel="Normalized power", xticks=[0, 5.0, 10.0, 15.0, 20.0])
        ax_right_tau.plot(np.arange(len(corr_right))*dt*1000, corr_right)
        plt.ylim([0, 0.1])
        plt.savefig(f"corr/right_corr_{os.path.basename(self.base_name)}.png")
        plt.close()
        fig = plt.figure(figsize=(14, 8))
        ax_left_tau = fig.add_subplot(
            111, xlabel=r'$\tau$ [ms]', ylabel="Normalized power", xticks=[0, 5.0, 10.0, 15.0, 20.0])
        ax_left_tau.plot(np.arange(len(corr_left))*dt*1000, corr_left)
        plt.ylim([0, 0.1])
        plt.savefig(f"corr/left_corr_{os.path.basename(self.base_name)}.png")
        plt.close()
        with open("./corr_{}/corr_{}.csv".format(os.path.basename(os.path.dirname(self.base_name)), os.path.splitext(self.base_name)[0].split("\\")[-1]), "a") as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(["tim", "echo_right", "echo_left",
                             "corr_right", "corr_left"])
            writer.writerows(data_for_csv)


def main(csv_list, emit_csv_list):
    """
    main for detect point from echo
    """
    for idx, csv in enumerate(csv_list):
        data_field = DataField(csv, emit_csv_list)
        print("analyze target:{}".format(csv))
        data_field.preprocessing()
        data_field.get_info()
        data_field.save()
        print("-------finish:{}th/{}-----------".format(idx+1, len(csv_list)))


if __name__ == "__main__":
    argvs = sys.argv
    if len(argvs) < 2:
        print(
            "Usage: python {} [folder of csv] [folder of emit_pulse csv]".format(argvs[0]))
        exit()
    csv_list = sorted(glob.glob("{}/*.csv".format(argvs[1])))
    if csv_list != []:
        if not os.path.exists("./corr_{}".format(os.path.basename(argvs[1]))):
            os.makedirs("./corr_{}".format(os.path.basename(argvs[1])))
        if not os.path.exists("./echo_point_{}".format(os.path.basename(argvs[1]))):
            os.makedirs("./echo_point_{}".format(os.path.basename(argvs[1])))
    print("csv_list:{}".format(csv_list))
    emit_csv_list = sorted(glob.glob("{}/*.csv".format(argvs[2])))
    main(csv_list, emit_csv_list)
