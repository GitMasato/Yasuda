from . import functions as func
import re
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib as path
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

temp_cmap = plt.get_cmap("Blues")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams["mathtext.fontset"] = 'dejavusans'
plt.rcParams['font.size'] = 10
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
cr_path = path.Path()


class Plural_data_output:

    def __init__(self):

        self.tgt_dirs = []
        self.output_name  = []
        self.ptcl_num = []
        self.ptcl_dia = []
        self.ptcl_coe_res = []
        self.end_excitation = []
        self.cell_length = []
        self.turn_on_for_accum = []
        self.turn_off_for_accum = []
        self.measurement = []
        self.normal_bins_size = []
        self.time_start_accum = []
        self.cubic_cell_bool = []
        self.vibration_on = []
        self.angle_xy = []
        self.angle_z = []
        self.vib_amplitude = []
        self.vib_frequency = []
        self.cmap_val = []
        self.ls_val = []
        self.lw_val = []
        self.label_val = []


    def set_values(self, Tgt_dirs, Output_name, Ptcl_num, Ptcl_dia, Ptcl_coe_res, End_excitation, Cell_length, Turn_on_for_accum, Turn_off_for_accum, Measurement, Normal_bins_size, Time_start_accum, Cubic_cell_bool, Vibration_on, Angle_xy, Angle_z, Vib_amplitude, Vib_frequency, Cmap_array, Ls_array, Lw_array, Label_array):

        for idx_dir in range (len(Tgt_dirs)):

            self.tgt_dirs.append(Tgt_dirs[idx_dir])

            if(len(Tgt_dirs) == len (Output_name)):
                self.output_name.append(Output_name[idx_dir])
            else:
                self.output_name.append(Output_name[0])

            if(len(Tgt_dirs) == len (Ptcl_num)):
                self.ptcl_num.append(Ptcl_num[idx_dir])
            else:
                self.ptcl_num.append(Ptcl_num[0])

            if(len(Tgt_dirs) == len (Ptcl_dia)):
                self.ptcl_dia.append(Ptcl_dia[idx_dir])
            else:
                self.ptcl_dia.append(Ptcl_dia[0])

            if(len(Tgt_dirs) == len (Ptcl_coe_res)):
                self.ptcl_coe_res.append(Ptcl_coe_res[idx_dir])
            else:
                self.ptcl_coe_res.append(Ptcl_coe_res[0])

            if(len(Tgt_dirs) == len (End_excitation)):
                self.end_excitation.append(End_excitation[idx_dir])
            else:
                self.end_excitation.append(End_excitation[0])

            if(len(Tgt_dirs) == len (Cell_length)):
                self.cell_length.append(Cell_length[idx_dir])
            else:
                self.cell_length.append(Cell_length[0])

            if(len(Tgt_dirs) == len (Turn_on_for_accum)):
                self.turn_on_for_accum.append(Turn_on_for_accum[idx_dir])
            else:
                self.turn_on_for_accum.append(Turn_on_for_accum[0])

            if(len(Tgt_dirs) == len (Turn_off_for_accum)):
                self.turn_off_for_accum.append(Turn_off_for_accum[idx_dir])
            else:
                self.turn_off_for_accum.append(Turn_off_for_accum[0])

            if(len(Tgt_dirs) == len (Measurement)):
                self.measurement.append(Measurement[idx_dir])
            else:
                self.measurement.append(Measurement[0])

            if(len(Tgt_dirs) == len (Normal_bins_size)):
                self.normal_bins_size.append(Normal_bins_size[idx_dir])
            else:
                self.normal_bins_size.append(Normal_bins_size[0])

            if(len(Tgt_dirs) == len (Time_start_accum)):
                self.time_start_accum.append(Time_start_accum[idx_dir])
            else:
                self.time_start_accum.append(Time_start_accum[0])

            if(len(Tgt_dirs) == len (Cubic_cell_bool)):
                self.cubic_cell_bool.append(Cubic_cell_bool[idx_dir])
            else:
                self.cubic_cell_bool.append(Cubic_cell_bool[0])

            if(len(Tgt_dirs) == len (Vibration_on)):
                self.vibration_on.append(Vibration_on[idx_dir])
            else:
                self.vibration_on.append(Vibration_on[0])

            if(len(Tgt_dirs) == len (Angle_xy)):
                self.angle_xy.append(Angle_xy[idx_dir])
            else:
                self.angle_xy.append(Angle_xy[0])

            if(len(Tgt_dirs) == len (Angle_z)):
                self.angle_z.append(Angle_z[idx_dir])
            else:
                self.angle_z.append(Angle_z[0])

            if(len(Tgt_dirs) == len (Vib_amplitude)):
                self.vib_amplitude.append(Vib_amplitude[idx_dir])
            else:
                self.vib_amplitude.append(Vib_amplitude[0])

            if(len(Tgt_dirs) == len (Vib_frequency)):
                self.vib_frequency.append(Vib_frequency[idx_dir])
            else:
                self.vib_frequency.append(Vib_frequency[0])

            if(len(Tgt_dirs) == len (Cmap_array)):
                self.cmap_val.append(Cmap_array[idx_dir])
            else:
                self.cmap_val.append(Cmap_array[0])

            if(len(Tgt_dirs) == len (Ls_array)):
                self.ls_val.append(Ls_array[idx_dir])
            else:
                self.ls_val.append(Ls_array[0])

            if(len(Tgt_dirs) == len (Lw_array)):
                self.lw_val.append(Lw_array[idx_dir])
            else:
                self.lw_val.append(Lw_array[0])

            if(len(Tgt_dirs) == len (Label_array)):
                self.label_val.append(Label_array[idx_dir])
            else:
                self.label_val.append(Label_array[0])


    def plural_ave_velo(self):

        plt.figure(figsize=(4.5,3),dpi=300)
        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            ave_velo_3D = np.zeros(len(csv_files))
            ave_velo_time = np.zeros(len(csv_files))

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                # read csv file
                df = pd.read_csv(val_csv)
                ave_velo_3D[idx_csv] = np.mean(np.sqrt((df['velo_x (m/s)'] * df['velo_x (m/s)']) + (df['velo_y (m/s)'] * df['velo_y (m/s)']) + (df['velo_z (m/s)'] * df['velo_z (m/s)'])))
                file_time = re.findall(r'\d+', str(val_csv.name))[0]
                ave_velo_time[idx_csv] = int(file_time)

            # plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw=self.lw_val[idx], label=self.label_val[idx])

            if idx == 0:

              plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c="b", ls=':', lw=1.6, label=self.label_val[idx])

            if idx == 1:

              plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c="g", ls='-', lw=1.2, label=self.label_val[idx])

            if idx == 2:

              plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c="r", ls='--', lw=1.4, label=self.label_val[idx])

            if idx == 3:

              plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c="c", ls='-.', lw=1.2, label=self.label_val[idx])

            if idx == 4:

              plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c="m", ls=':', lw=1.2, label=self.label_val[idx])

        plt.xlabel('time  $\it{s}$')
        plt.ylabel('average velocity  $\it{mm/s}$')
        # plt.xlim(xmin=0, xmax=3.0)
        # plt.xticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], [0, 0.5, 1.0, 1.5, 2.0, 2.5,
        # 3.0], rotation=0)
        plt.xlim(xmin=0, xmax=5)
        plt.xticks([0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], rotation=0)

        # plt.ylim(ymin=0,ymax=50)
        # plt.yticks([0, 10, 20, 30, 40, 50], [0, 10, 20, 30, 40, 50],rotation=0)
        # plt.ylim(ymin=0,ymax=30)
        # plt.yticks([0, 5, 10, 15, 20, 25, 30], [0, 5, 10, 15, 20, 25, 30],rotation=0)
        # plt.ylim(ymin=0,ymax=70)
        # plt.yticks([0, 10, 20, 30, 40, 50, 60, 70], [0, 10, 20, 30, 40, 50, 60,
        # 70],rotation=0)
        plt.ylim(ymin=0,ymax=140)
        plt.yticks([0, 20, 40, 60, 80, 100, 120, 140], [0, 20, 40, 60, 80, 100, 120, 140],rotation=0)


        if(len(self.tgt_dirs) >= 3):
            # plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0, ncol=2)
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0, ncol=2)
        else:
            plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0)
            # plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_plural_ave_velo.png'.format(str(self.output_name[idx])),bbox_inches='tight')
        #plt.show()
        plt.close()


    def plural_saturation_velo(self):

        with open('_post/_plural_data/{0}_plural_saturation_velo.csv'.format(str(self.output_name[0])), 'w') as f:

            writer = csv.writer(f)
            writer.writerow(['name', 'turn_on_for_accum_ms', 'ref_time_ms', 'ref_ave_velo_mm/s', 'time_ms', 'ave_velo_mm/s'])

            # search and process in each directory
            for idx in range(len(self.tgt_dirs)):

                out_put_list = []
                csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
                csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))

                df_ref = pd.read_csv(csv_files_path / '_{0:0>5}_ms.csv'.format(int(self.end_excitation[idx]*1000)))
                ref_velo = np.mean(np.sqrt((df_ref['velo_x (m/s)'] * df_ref['velo_x (m/s)']) + (df_ref['velo_y (m/s)'] * df_ref['velo_y (m/s)']) + (df_ref['velo_z (m/s)'] * df_ref['velo_z (m/s)'])))
                out_put_list.append(self.label_val[idx])
                out_put_list.append(self.turn_on_for_accum[idx]*1000)
                out_put_list.append(self.end_excitation[idx]*1000)
                out_put_list.append(ref_velo*1000)
                temp_flag = True

                # search and process each csv file
                for idx_csv, val_csv in enumerate (csv_files):

                    if(temp_flag==True):

                        # read csv file
                        df = pd.read_csv(val_csv)
                        ave_velo = np.mean(np.sqrt((df['velo_x (m/s)'] * df['velo_x (m/s)']) + (df['velo_y (m/s)'] * df['velo_y (m/s)']) + (df['velo_z (m/s)'] * df['velo_z (m/s)'])))

                        if(ave_velo >= ref_velo):
                            out_put_list.append(int(re.findall(r'\d+', str(val_csv.name))[0]))
                            out_put_list.append(ave_velo*1000)
                            writer.writerow(out_put_list)
                            temp_flag = False

        df_csv = pd.read_csv('_post/_plural_data/{0}_plural_saturation_velo.csv'.format(str(self.output_name[0])), header=0)
        df_1 = df_csv[df_csv['name'].str.contains('0.5x')]
        df_2 = df_csv[df_csv['name'].str.contains('1x')]
        df_3 = df_csv[df_csv['name'].str.contains('2x')]
        temp_x_1 = np.linspace(np.min(df_1['turn_on_for_accum_ms']), np.max(df_1['turn_on_for_accum_ms']))
        temp_x_2 = np.linspace(np.min(df_2['turn_on_for_accum_ms']), np.max(df_2['turn_on_for_accum_ms']))
        temp_x_3 = np.linspace(np.min(df_3['turn_on_for_accum_ms']), np.max(df_3['turn_on_for_accum_ms']))


        plt.figure(figsize=(4.5,3),dpi=300)
        plt.plot(df_1['turn_on_for_accum_ms'], df_1['time_ms'] * 0.001, 's', color='none', markersize=5, markeredgewidth=2, markeredgecolor='blue', alpha=0.8, label='0.5x')
        plt.plot(df_2['turn_on_for_accum_ms'], df_2['time_ms'] * 0.001, '^', color='none', markersize=5, markeredgewidth=2, markeredgecolor='green', alpha=0.8, label='1x')
        plt.plot(df_3['turn_on_for_accum_ms'], df_3['time_ms'] * 0.001, 'o', color='none', markersize=5, markeredgewidth=2, markeredgecolor='red', alpha=0.8, label='2x')
        param_1, cov_1 = curve_fit(func.exp_decay_fit, df_1['turn_on_for_accum_ms'], df_1['time_ms'] * 0.001)
        param_2, cov_2 = curve_fit(func.exp_decay_fit, df_2['turn_on_for_accum_ms'], df_2['time_ms'] * 0.001)
        param_3, cov_3 = curve_fit(func.exp_decay_fit, df_3['turn_on_for_accum_ms'], df_3['time_ms'] * 0.001)
        plt.plot(temp_x_1, func.exp_decay_fit(temp_x_1, *param_1), color='blue', ls=':', lw = 1.6)
        plt.plot(temp_x_2, func.exp_decay_fit(temp_x_2, *param_2), color='green', ls='-', lw = 1.2)
        plt.plot(temp_x_3, func.exp_decay_fit(temp_x_3, *param_3), color='red', ls='--', lw = 1.4)
        plt.xlabel('turn-on duration  $\it{ms}$')
        plt.ylabel('saturation time  $\it{s}$')
        #plt.xlim(xmin=0, xmax=100)
        #plt.ylim(ymin=0,ymax=2.2)
        #plt.xticks([0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], rotation=0)
        #plt.yticks([0, 25, 50, 75], [0, 25, 50, 75],rotation=0)
        plt.legend(bbox_to_anchor=(0.8, 0.5), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_plural_saturation_time.png'.format(str(self.output_name[0])), bbox_inches='tight')
        #plt.show()
        plt.close()


        plt.figure(figsize=(4.5,3),dpi=300)
        plt.plot(df_1['turn_on_for_accum_ms'], df_1['ref_ave_velo_mm/s'], 's', color='none', markersize=3, markeredgewidth=2, markeredgecolor='blue', alpha=0.8, label='0.5x')
        plt.plot(df_2['turn_on_for_accum_ms'], df_2['ref_ave_velo_mm/s'], '^', color='none', markersize=3, markeredgewidth=2, markeredgecolor='green', alpha=0.8, label='1x')
        plt.plot(df_3['turn_on_for_accum_ms'], df_3['ref_ave_velo_mm/s'], 'o', color='none', markersize=3, markeredgewidth=2, markeredgecolor='red', alpha=0.8, label='2x')

        # param_1, cov_1 = curve_fit(func.log_fit, df_1['turn_on_for_accum_ms'], df_1['ref_ave_velo_mm/s'], bounds=((1.0, 0.1, -10.0),(100.0, 100.0, 100.0)))
        # param_2, cov_2 = curve_fit(func.log_fit, df_2['turn_on_for_accum_ms'], df_2['ref_ave_velo_mm/s'], bounds=((1.0, 0.1, -10.0),(100.0, 100.0, 100.0)))
        # param_3, cov_3 = curve_fit(func.log_fit, df_3['turn_on_for_accum_ms'], df_3['ref_ave_velo_mm/s'], bounds=((1.0, 0.1, -10.0),(100.0, 100.0, 100.0)))

        # param_1, cov_1 = curve_fit(lambda t,a,b: a+b*np.log(t), df_1['turn_on_for_accum_ms'][:8], df_1['ref_ave_velo_mm/s'][:8])
        # param_2, cov_2 = curve_fit(lambda t,a,b: a+b*np.log(t), df_2['turn_on_for_accum_ms'][:7], df_2['ref_ave_velo_mm/s'][:7])
        # param_3, cov_3 = curve_fit(lambda t,a,b: a+b*np.log(t), df_3['turn_on_for_accum_ms'][:6], df_3['ref_ave_velo_mm/s'][:6])
        # plt.plot(temp_x_1, param_1[0] + param_1[1] * np.log(temp_x_1), color='blue', ls=':', lw = 1.0)
        # plt.plot(temp_x_2, param_2[0] + param_2[1] * np.log(temp_x_1), color='green', ls=':', lw = 1.0)
        # plt.plot(temp_x_3, param_3[0] + param_3[1] * np.log(temp_x_1), color='red', ls=':', lw = 1.0)

        param_1, cov_1 = curve_fit(lambda t,a,b: a+b*t, df_1['turn_on_for_accum_ms'][1:8], df_1['ref_ave_velo_mm/s'][1:8])
        param_2, cov_2 = curve_fit(lambda t,a,b: a+b*t, df_2['turn_on_for_accum_ms'][1:8], df_2['ref_ave_velo_mm/s'][1:8])
        param_3, cov_3 = curve_fit(lambda t,a,b: a+b*t, df_3['turn_on_for_accum_ms'][1:7], df_3['ref_ave_velo_mm/s'][1:7])
        plt.plot(temp_x_1, (param_1[0] - 0) + (param_1[1] + 0.0) * temp_x_1, color='blue', ls=':', lw = 1.0)
        plt.plot(temp_x_2, (param_2[0] - 0) + (param_2[1] + 0.02) * temp_x_1, color='green', ls=':', lw = 1.0)
        plt.plot(temp_x_3, (param_3[0] - 1) + (param_3[1] + 0.0) * temp_x_1, color='red', ls=':', lw = 1.0)

        plt.xlabel('turning-on duration  $\it{ms}$')
        plt.ylabel('average velocity  $\it{mm/s}$')
        plt.xlim(xmin=-1, xmax=102)
        plt.xticks([0, 20, 40, 60, 80, 100], [0, 20, 40, 60, 80, 100], rotation=0)
        plt.ylim(ymin=0,ymax=140)
        plt.yticks([0, 20, 40, 60, 80, 100, 120, 140], [0, 20, 40, 60, 80, 100, 120, 140],rotation=0)

        # plt.text(80, 10, "{0:.2f}+{1:.2f}".format(param_1[0], param_1[1]) + 'log($\it{t}_{E}$)', va="center", ha="center", color='blue', rotation=0)
        # plt.text(80, 35, "{0:.2f}+{1:.2f}".format(param_2[0], param_2[1]) + 'log($\it{t}_{E}$)', va="center", ha="center", color='green', rotation=0)
        # plt.text(30, 70, "{0:.2f}+{1:.2f}".format(param_3[0], param_3[1]) + 'log($\it{t}_{E}$)', va="center", ha="center", color='red', rotation=0)

        plt.legend(bbox_to_anchor=(0.3, 0.6), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_plural_saturation_velo.png'.format(str(self.output_name[0])), bbox_inches='tight')
        #plt.show()
        plt.close()



    def plural_cell_layer_velo_dist(self):

        plt.figure(figsize=(4.5,3),dpi=300)
        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
            half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
            velo_dist_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), velo_dist_measurement))
            df_pos = pd.DataFrame()

            if(self.cubic_cell_bool[idx]):
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
            else:
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03

            if((self.vibration_on[idx]) and (velo_dist_measurement <= int((self.end_excitation[idx] * 1000) + 0.5))):

                temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(velo_dist_measurement) * 0.001)
                df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(velo_dist_measurement) * 0.001)
                df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(velo_dist_measurement) * 0.001)


            if(self.cubic_cell_bool[idx]):

                df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                ave_velo_in_each_cell = np.zeros(half_cubic_length)
                ptcl_num_in_each_cell = np.zeros(half_cubic_length)

                for idx_ptcl in range (len(df_pos)):

                    temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                    dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<=half_cubic_length - 1) else (temp_pos[0] - half_cubic_length).astype(int)
                    dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<=half_cubic_length - 1) else (temp_pos[1] - half_cubic_length).astype(int)
                    dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<=half_cubic_length - 1) else (temp_pos[2] - half_cubic_length).astype(int)
                    in_cell = max([dis_x, dis_y, dis_z])
                    ave_velo_in_each_cell[in_cell] += df.at[idx_ptcl, 'velo_ab (m/s)']
                    ptcl_num_in_each_cell[in_cell] += 1

                for idx_cell in range (4):
                    ave_velo_in_each_cell[4] += ave_velo_in_each_cell[idx_cell]
                    ptcl_num_in_each_cell[4] += ptcl_num_in_each_cell[idx_cell]
                    ave_velo_in_each_cell[idx_cell] = 0
                    ptcl_num_in_each_cell[idx_cell] = 0

                output_ave_velo_in_each_cell = 1000 * np.divide(ave_velo_in_each_cell, ptcl_num_in_each_cell, out=np.zeros_like(ave_velo_in_each_cell), where=ptcl_num_in_each_cell!=0)
                plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_velo_in_each_cell, c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw=self.lw_val[idx], label=self.label_val[idx])

            else:
                rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)
                ave_velo_in_each_cell = np.zeros(sphere_radius)
                ptcl_num_in_each_cell = np.zeros(sphere_radius)

                for idx_ptcl in range (len(df_pos)):
                    ave_velo_in_each_cell[rel_pos[idx_ptcl]] += df.at[idx_ptcl, 'velo_ab (m/s)']
                    ptcl_num_in_each_cell[rel_pos[idx_ptcl]] += 1

                for idx_cell in range (4):
                    ave_velo_in_each_cell[4] += ave_velo_in_each_cell[idx_cell]
                    ptcl_num_in_each_cell[4] += ptcl_num_in_each_cell[idx_cell]
                    ave_velo_in_each_cell[idx_cell] = 0
                    ptcl_num_in_each_cell[idx_cell] = 0

                output_ave_velo_in_each_cell = 1000 * np.divide(ave_velo_in_each_cell, ptcl_num_in_each_cell, out=np.zeros_like(ave_velo_in_each_cell), where=ptcl_num_in_each_cell!=0)
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw=self.lw_val[idx], label=self.label_val[idx])

        plt.xlabel('layer')
        plt.ylabel('velocity  $\it{mm/s}$')
        plt.xlim(xmin=5, xmax=25)
        plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
        #plt.ylim(ymin=0,ymax=50.0)
        #plt.yticks([0, 10, 20, 30, 40, 50], [0, 10, 20, 30, 40, 50],rotation=0)

        if(len(self.tgt_dirs) >= 3):
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0, ncol=2)
        else:
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_plural_velo_cell_{1}_ms.png'.format(str(self.output_name[idx]), velo_dist_measurement), bbox_inches='tight')
        #plt.show()
        plt.close()


    def plural_accumulated_cell_layer_velo_dist(self):

        plt.figure(figsize=(4.5,3),dpi=300)
        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count = 0

            if(self.cubic_cell_bool[idx]):
                half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
                ave_velo_in_each_cell = np.zeros(half_cubic_length)
                ptcl_num_in_each_cell = np.zeros(half_cubic_length)
            else:
                sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
                ave_velo_in_each_cell = np.zeros(sphere_radius)
                ptcl_num_in_each_cell = np.zeros(sphere_radius)

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):

                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        df_pos = pd.DataFrame()
                        count += 1

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
                        else:
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03   # adjust for calculations condition

                        if((self.vibration_on[idx]) and (int(file_time) <= int(self.end_excitation[idx] * 1000))):

                            temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                            vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                            df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                        else:
                            rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)


                        for idx_ptcl in range (len(df_pos)):
                            if(self.cubic_cell_bool[idx]):
                                temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                                dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<=half_cubic_length - 1) else (temp_pos[0] - half_cubic_length).astype(int)
                                dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<=half_cubic_length - 1) else (temp_pos[1] - half_cubic_length).astype(int)
                                dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<=half_cubic_length - 1) else (temp_pos[2] - half_cubic_length).astype(int)
                                in_cell = max([dis_x, dis_y, dis_z])
                                ave_velo_in_each_cell[in_cell] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                ptcl_num_in_each_cell[in_cell] += 1
                            else:
                                ave_velo_in_each_cell[rel_pos[idx_ptcl]] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                ptcl_num_in_each_cell[rel_pos[idx_ptcl]] += 1

            for idx_cell in range (4):
                ave_velo_in_each_cell[4] += ave_velo_in_each_cell[idx_cell]
                ptcl_num_in_each_cell[4] += ptcl_num_in_each_cell[idx_cell]
                ave_velo_in_each_cell[idx_cell] = 0
                ptcl_num_in_each_cell[idx_cell] = 0
            output_ave_velo_in_each_cell = 1000 * np.divide(ave_velo_in_each_cell, ptcl_num_in_each_cell, out=np.zeros_like(ave_velo_in_each_cell), where=ptcl_num_in_each_cell!=0)

            if(self.cubic_cell_bool[idx]):
                # plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_velo_in_each_cell, c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw=self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="m", ls=':', lw=1.2, label=self.label_val[idx])

            else:
                # plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw=self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="m", ls=':', lw=1.2, label=self.label_val[idx])

        plt.xlabel('layer')
        plt.ylabel('velocity  $\it{mm/s}$')
        plt.xlim(xmin=5, xmax=25)
        #plt.ylim(ymin=0,ymax=140.0)
        plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
        #plt.ylim(ymin=0,ymax=50.0)
        #plt.yticks([0, 10, 20, 30, 40, 50], [0, 10, 20, 30, 40, 50],rotation=0)

        if(len(self.tgt_dirs) >= 3):
            # plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0, ncol=2)
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0, ncol=2)
        else:
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_{1}_accum_{2}_ptcls_end_excitation_plural_velo_cell.png'.format(str(self.output_name[idx]), count, count * len(df)), bbox_inches='tight')
        #plt.show()
        plt.close()


        plt.figure(figsize=(4.5,3),dpi=300)
        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count = 0

            if(self.cubic_cell_bool[idx]):
                half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
                ave_velo_in_each_cell = np.zeros(half_cubic_length)
                ptcl_num_in_each_cell = np.zeros(half_cubic_length)
            else:
                sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
                ave_velo_in_each_cell = np.zeros(sphere_radius)
                ptcl_num_in_each_cell = np.zeros(sphere_radius)

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):

                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        df_pos = pd.DataFrame()
                        count += 1

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
                        else:
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03   # adjust for calculations condition

                        if((self.vibration_on[idx]) and (int(file_time) <= int(self.end_excitation[idx] * 1000))):

                            temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                            vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                            df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                        else:
                            rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)


                        for idx_ptcl in range (len(df_pos)):
                            if(self.cubic_cell_bool[idx]):
                                temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                                dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<=half_cubic_length - 1) else (temp_pos[0] - half_cubic_length).astype(int)
                                dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<=half_cubic_length - 1) else (temp_pos[1] - half_cubic_length).astype(int)
                                dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<=half_cubic_length - 1) else (temp_pos[2] - half_cubic_length).astype(int)
                                in_cell = max([dis_x, dis_y, dis_z])
                                ave_velo_in_each_cell[in_cell] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                ptcl_num_in_each_cell[in_cell] += 1
                            else:
                                ave_velo_in_each_cell[rel_pos[idx_ptcl]] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                ptcl_num_in_each_cell[rel_pos[idx_ptcl]] += 1

            for idx_cell in range (4):
                ave_velo_in_each_cell[4] += ave_velo_in_each_cell[idx_cell]
                ptcl_num_in_each_cell[4] += ptcl_num_in_each_cell[idx_cell]
                ave_velo_in_each_cell[idx_cell] = 0
                ptcl_num_in_each_cell[idx_cell] = 0
            output_ave_velo_in_each_cell = 1000 * np.divide(ave_velo_in_each_cell, ptcl_num_in_each_cell, out=np.zeros_like(ave_velo_in_each_cell), where=ptcl_num_in_each_cell!=0)

            if(self.cubic_cell_bool[idx]):
                # plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_velo_in_each_cell, c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw=self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c="m", ls=':', lw=1.2, label=self.label_val[idx])

            else:
                # plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw=self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c="m", ls=':', lw=1.2, label=self.label_val[idx])

        plt.xlabel('layer')
        plt.ylabel('velocity  $\it{mm/s}$')
        plt.xlim(xmin=5, xmax=25)
        plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)

        # plt.ylim(ymin=0,ymax=50)
        # plt.yticks([0, 10, 20, 30, 40, 50], [0, 10, 20, 30, 40, 50],rotation=0)
        plt.ylim(ymin=0,ymax=30)
        plt.yticks([0, 5, 10, 15, 20, 25, 30], [0, 5, 10, 15, 20, 25, 30],rotation=0)
        # plt.ylim(ymin=0,ymax=70)
        # plt.yticks([0, 10, 20, 30, 40, 50, 60, 70], [0, 10, 20, 30, 40, 50, 60,
        # 70],rotation=0)
        # plt.ylim(ymin=0,ymax=140)
        # plt.yticks([0, 20, 40, 60, 80, 100, 120, 140], [0, 20, 40, 60, 80, 100, 120, 140],rotation=0)

        if(len(self.tgt_dirs) >= 3):
            # plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0, ncol=2)
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0, ncol=2)
        else:
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_{1}_accum_{2}_ptcls_end_cooling_plural_velo_cell.png'.format(str(self.output_name[idx]), count, count * len(df)), bbox_inches='tight')
        #plt.show()
        plt.close()


    def plural_cell_layer_pos_dist(self):

        plt.figure(figsize=(4.5,3),dpi=300)
        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            pos_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), pos_measurement))
            df_pos = pd.DataFrame()
            half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
            sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)

            if(self.cubic_cell_bool[idx]):
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition

            else:
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03

            if((self.vibration_on[idx]) and (pos_measurement <= int((self.end_excitation[idx] * 1000) + 0.5))):

                temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(pos_measurement) * 0.001)
                df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(pos_measurement) * 0.001)
                df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(pos_measurement) * 0.001)


            if(self.cubic_cell_bool[idx]):

                df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                cell_number = np.append(np.array([2**3]), np.diff(np.power(np.arange(2, 50+2, 2), 3)))
                ptcl_num_in_each_cell = np.zeros(half_cubic_length)

                # average velocity in each cell layer at each tine step
                for idx_ptcl in range (len(df_pos)):

                    temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                    dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<half_cubic_length) else (temp_pos[0] - half_cubic_length).astype(int)
                    dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<half_cubic_length) else (temp_pos[1] - half_cubic_length).astype(int)
                    dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<half_cubic_length) else (temp_pos[2] - half_cubic_length).astype(int)
                    in_cell = max([dis_x, dis_y, dis_z])
                    ptcl_num_in_each_cell[in_cell] += 1

                for idx_cell in range (4):
                    ptcl_num_in_each_cell[4] += ptcl_num_in_each_cell[idx_cell]
                    cell_number[4] += cell_number[idx_cell]
                    ptcl_num_in_each_cell[idx_cell] = 0
                    cell_number[idx_cell] = 0

                output_ave_pos_in_each_cell = np.divide(ptcl_num_in_each_cell, cell_number, out=np.zeros_like(ptcl_num_in_each_cell), where=ptcl_num_in_each_cell!=0)
                plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw = self.lw_val[idx], label=self.label_val[idx])

            else:
                rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)
                ptcl_num_in_each_layer = np.zeros(sphere_radius)
                volume_in_each_layer = np.diff(4.0 / 3.0 * np.pi * np.power(np.array(range(0,sphere_radius+1)),3))

                for idx_ptcl in range (len(df_pos)):
                    ptcl_num_in_each_layer[rel_pos[idx_ptcl]] += 1

                for idx_cell in range (4):
                    ptcl_num_in_each_layer[4] += ptcl_num_in_each_layer[idx_cell]
                    volume_in_each_layer[4] += volume_in_each_layer[idx_cell]
                    ptcl_num_in_each_layer[idx_cell] = 0
                    volume_in_each_layer[idx_cell] = 0

                output_ave_pos_in_each_layer = np.divide(ptcl_num_in_each_layer, volume_in_each_layer, out=np.zeros_like(ptcl_num_in_each_layer), where=ptcl_num_in_each_layer!=0)
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_layer/np.sum(output_ave_pos_in_each_layer), c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw = self.lw_val[idx], label=self.label_val[idx])

        plt.xlabel('layer')
        plt.ylabel('probability density')
        plt.xlim(xmin=5,xmax=25)
        plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
        plt.ylim(ymin=0,ymax=0.1)
        plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], [0, 0.02, 0.04, 0.06, 0.08, 0.1],rotation=0)

        if(len(self.tgt_dirs) >= 3):
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0, ncol=2)
        else:
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_plural_pos_cell_{1}_ms.png'.format(str(self.output_name[idx]), pos_measurement), bbox_inches='tight')
        #plt.show()
        plt.close()


    def plural_accumulated_cell_layer_pos_dist(self):

        plt.figure(figsize=(4.5,3),dpi=300)
        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count = 0
            half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
            cell_number = np.append(np.array([2**3]), np.diff(np.power(np.arange(2, 50+2, 2), 3)))
            sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
            volume_in_each_cell = np.diff(4.0 / 3.0 * np.pi * np.power(np.array(range(0,sphere_radius+1)),3))
            if(self.cubic_cell_bool[idx]):
                ptcl_num_in_each_cell = np.zeros(half_cubic_length)
            else:
                ptcl_num_in_each_cell = np.zeros(sphere_radius)


            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):

                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        df_pos = pd.DataFrame()
                        count += 1

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
                        else:
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03   # adjust for calculations condition

                        if((self.vibration_on[idx]) and (int(file_time) <= int(self.end_excitation[idx] * 1000))):

                            temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                            vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                            df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                        else:
                            rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)


                        for idx_ptcl in range (len(df_pos)):
                            if(self.cubic_cell_bool[idx]):
                                temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                                dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<=half_cubic_length - 1) else (temp_pos[0] - half_cubic_length).astype(int)
                                dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<=half_cubic_length - 1) else (temp_pos[1] - half_cubic_length).astype(int)
                                dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<=half_cubic_length - 1) else (temp_pos[2] - half_cubic_length).astype(int)
                                in_cell = max([dis_x, dis_y, dis_z])
                                ptcl_num_in_each_cell[in_cell] += 1
                            else:
                                ptcl_num_in_each_cell[rel_pos[idx_ptcl]] += 1


            for idx_cell in range (4):
                ptcl_num_in_each_cell[4] += ptcl_num_in_each_cell[idx_cell]
                cell_number[4] += cell_number[idx_cell]
                volume_in_each_cell[4] += volume_in_each_cell[idx_cell]
                ptcl_num_in_each_cell[idx_cell] = 0
                cell_number[idx_cell] = 0
                volume_in_each_cell[idx_cell] = 0

            if(self.cubic_cell_bool[idx]):
                output_ave_pos_in_each_cell = np.divide(ptcl_num_in_each_cell, cell_number, out=np.zeros_like(ptcl_num_in_each_cell), where=ptcl_num_in_each_cell!=0)

                # plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw = self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="m", ls=':', lw=1.2, label=self.label_val[idx])

            else:
                output_ave_pos_in_each_cell = np.divide(ptcl_num_in_each_cell, volume_in_each_cell, out=np.zeros_like(ptcl_num_in_each_cell), where=ptcl_num_in_each_cell!=0)

                # plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw = self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,sphere_radius + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,sphere_radius + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="m", ls=':', lw=1.2, label=self.label_val[idx])


        plt.xlabel('layer')
        plt.ylabel('probability')
        plt.xlim(xmin=5,xmax=25)
        plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
        plt.ylim(ymin=0,ymax=0.1)
        plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], [0, 0.02, 0.04, 0.06, 0.08, 0.1],rotation=0)

        if(len(self.tgt_dirs) >= 3):
            # plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0, ncol=2)
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0, ncol=2)
        else:
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0)
            # plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0)

        plt.savefig('_post/_plural_data/{0}_{1}_accum_{2}_ptcls_end_excitation_plural_pos_cell.png'.format(str(self.output_name[idx]), count, count * len(df)), bbox_inches='tight')
        #plt.show()
        plt.close()


        plt.figure(figsize=(4.5,3),dpi=300)
        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count = 0
            half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
            cell_number = np.append(np.array([2**3]), np.diff(np.power(np.arange(2, 50+2, 2), 3)))
            sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
            volume_in_each_cell = np.diff(4.0 / 3.0 * np.pi * np.power(np.array(range(0,sphere_radius+1)),3))
            if(self.cubic_cell_bool[idx]):
                ptcl_num_in_each_cell = np.zeros(half_cubic_length)
            else:
                ptcl_num_in_each_cell = np.zeros(sphere_radius)


            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):

                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        df_pos = pd.DataFrame()
                        count += 1

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
                        else:
                            df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03   # adjust for calculations condition

                        if((self.vibration_on[idx]) and (int(file_time) <= int(self.end_excitation[idx] * 1000))):

                            temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                            vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                            vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                            df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)

                        if(self.cubic_cell_bool[idx]):
                            df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                        else:
                            rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)


                        for idx_ptcl in range (len(df_pos)):
                            if(self.cubic_cell_bool[idx]):
                                temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                                dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<=half_cubic_length - 1) else (temp_pos[0] - half_cubic_length).astype(int)
                                dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<=half_cubic_length - 1) else (temp_pos[1] - half_cubic_length).astype(int)
                                dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<=half_cubic_length - 1) else (temp_pos[2] - half_cubic_length).astype(int)
                                in_cell = max([dis_x, dis_y, dis_z])
                                ptcl_num_in_each_cell[in_cell] += 1
                            else:
                                ptcl_num_in_each_cell[rel_pos[idx_ptcl]] += 1

            for idx_cell in range (4):
                ptcl_num_in_each_cell[4] += ptcl_num_in_each_cell[idx_cell]
                cell_number[4] += cell_number[idx_cell]
                volume_in_each_cell[4] += volume_in_each_cell[idx_cell]
                ptcl_num_in_each_cell[idx_cell] = 0
                cell_number[idx_cell] = 0
                volume_in_each_cell[idx_cell] = 0

            if(self.cubic_cell_bool[idx]):
                output_ave_pos_in_each_cell = np.divide(ptcl_num_in_each_cell, cell_number, out=np.zeros_like(ptcl_num_in_each_cell), where=ptcl_num_in_each_cell!=0)

                # plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw = self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="m", ls=':', lw=1.2, label=self.label_val[idx])

            else:
                output_ave_pos_in_each_cell = np.divide(ptcl_num_in_each_cell, volume_in_each_cell, out=np.zeros_like(ptcl_num_in_each_cell), where=ptcl_num_in_each_cell!=0)

                # plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c=temp_cmap(self.cmap_val[idx]), ls=self.ls_val[idx], lw = self.lw_val[idx], label=self.label_val[idx])

                if idx == 0:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="b", ls=':', lw=1.6, label=self.label_val[idx])

                if idx == 1:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="g", ls='-', lw=1.2, label=self.label_val[idx])

                if idx == 2:

                  plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="r", ls='--', lw=1.4, label=self.label_val[idx])

                if idx == 3:

                  plt.plot(np.array(range(1,sphere_radius + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="c", ls='-.', lw=1.2, label=self.label_val[idx])

                if idx == 4:

                  plt.plot(np.array(range(1,sphere_radius + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c="m", ls=':', lw=1.2, label=self.label_val[idx])

        plt.xlabel('layer')
        plt.ylabel('probability')
        plt.xlim(xmin=5,xmax=25)
        plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
        plt.ylim(ymin=0,ymax=0.1)
        plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], [0, 0.02, 0.04, 0.06, 0.08, 0.1],rotation=0)

        if(len(self.tgt_dirs) >= 3):
            # plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0, ncol=2)
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0, ncol=2)
        else:
            # plt.legend(bbox_to_anchor=(0.8, 0.1), loc='lower right', borderaxespad=0)
            plt.legend(bbox_to_anchor=(1.0, 1.01), loc='lower right', borderaxespad=0)
        plt.savefig('_post/_plural_data/{0}_{1}_accum_{2}_ptcls_end_cooling_plural_pos_cell.png'.format(str(self.output_name[idx]), count, count * len(df)), bbox_inches='tight')
        #plt.show()
        plt.close()

