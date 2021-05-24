import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib as path
import scipy.stats as sp
from scipy.stats import norm
from scipy.stats import maxwell

temp_cmap = plt.get_cmap("Blues")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams["mathtext.fontset"] = 'dejavusans'
plt.rcParams['font.size'] = 10
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
cr_path = path.Path()


class Velo_data_output:

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
                self. angle_z.append(Angle_z[0])

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


    def ave_velo(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            ave_velo_3D = np.zeros(len(csv_files)) 
            ave_velo_2D = np.zeros(len(csv_files))  
            ave_velo_time = np.zeros(len(csv_files))
            sphere_cell_radius = self.cell_length[idx] * 0.5

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                # read csv file
                df = pd.read_csv(val_csv)
                ave_velo_3D[idx_csv] = np.mean(np.sqrt((df['velo_x (m/s)'] * df['velo_x (m/s)']) + (df['velo_y (m/s)'] * df['velo_y (m/s)']) + (df['velo_z (m/s)'] * df['velo_z (m/s)'])))
                ave_velo_2D[idx_csv] = np.mean(np.sqrt((df['velo_y (m/s)'] * df['velo_y (m/s)']) + (df['velo_z (m/s)'] * df['velo_z (m/s)'])))
                file_time = re.findall(r'\d+', str(val_csv.name))[0]
                ave_velo_time[idx_csv] = int(file_time)        

            # showing particle average-velocity
            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c=temp_cmap(1.0), ls='solid', lw = 0.7, label='calculation (3D)')
            plt.plot(ave_velo_time / 1000, ave_velo_2D * 1000, c=temp_cmap(0.5), ls='solid', lw = 0.7, label='calculation (2D)')
            plt.xlabel('time  $\it{s}$')
            plt.ylabel('average velocity  $\it{mm/s}$')        
            plt.xlim(xmin=0, xmax=self.end_excitation[idx])
            #plt.xticks([0, 15, 30, 45], [0, 15, 30, 45], rotation=0)
            #plt.ylim(ymin=0,ymax=50)
            #plt.yticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)
            #plt.grid(which='minor')
            plt.legend(bbox_to_anchor=(1, 1.01), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_ave_velo_1.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx])),bbox_inches='tight')
            #plt.show()
            plt.close()


            # showing particle average-velocity
            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(ave_velo_time / 1000, ave_velo_3D * 1000, c=temp_cmap(1.0), ls='solid', lw = 0.7, label='calculation (3D)')
            plt.plot(ave_velo_time / 1000, ave_velo_2D * 1000, c=temp_cmap(0.5), ls='solid', lw = 0.7, label='calculation (2D)')

            # showing haff's law   
            if(self.cubic_cell_bool[idx]):
                ptcl_num_den = self.ptcl_num[idx] / pow((self.cell_length[idx]*1000),3)            
            else:
                ptcl_num_den = self.ptcl_num[idx] / ((4/3) * np.pi * pow((sphere_cell_radius*1000),3))
            
            haff_ini = float(ave_velo_3D[(np.where(ave_velo_time / 1000 == self.end_excitation[idx]))[0]]) * 1000
            tau_inv = np.pi * ptcl_num_den * pow((self.ptcl_dia[idx]*1000),2) * (1 - pow(self.ptcl_coe_res[idx], 2)) * haff_ini * 0.5
            haff_x = np.arange(self.end_excitation[idx], self.end_excitation[idx] + 50.0, 0.01)
            haff_y = haff_ini / (1.0 + ((haff_x - self.end_excitation[idx]) * tau_inv))
            plt.plot(haff_x, haff_y, label="$\it{Haff}$'s theory (3D, $\it{\u03c4}$"+'= {:.2f})'.format(1.0/tau_inv), c='r', ls='dashed', lw = 0.7)
            plt.xlabel('time  $\it{s}$')
            plt.ylabel('average velocity  $\it{mm/s}$')        
            plt.xlim(xmin=int(self.end_excitation[idx]+0.5) - 1, xmax=int(self.end_excitation[idx]+0.5) + 5)
            #plt.xticks([4, 6, 8, 10], [4, 6, 8, 10],rotation=0)
            #plt.ylim(ymin=0,ymax=50)
            #plt.yticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)
            #plt.gca().yaxis.set_minor_locator(tick.MultipleLocator(10))
            #plt.grid(which='minor')
            plt.legend(bbox_to_anchor=(1, 1.01), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_ave_velo_2.png'.format(str(self.tgt_dirs[idx]),str(self.output_name[idx])),bbox_inches='tight')
            #plt.show()
            plt.close()
        

    def velo_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            velo_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), velo_measurement))        
            velo_x = df['velo_x (m/s)']
            velo_y = df['velo_y (m/s)']
            velo_z = df['velo_z (m/s)']

            std_velo_dis_x = (velo_x - np.mean(velo_x)) / np.std(velo_x)
            #moment1_x = np.mean(std_velo_dis_x)
            #moment2_x = np.var(std_velo_dis_x)
            moment3_x = sp.skew(std_velo_dis_x)
            moment4_x = sp.kurtosis(std_velo_dis_x)
            hist_x = plt.hist(std_velo_dis_x, bins=self.normal_bins_size[idx], density=True, log=True)
            plt.close()

            std_velo_dis_y = (velo_y - np.mean(velo_y)) / np.std(velo_y)
            #moment1_y = np.mean(std_velo_dis_y)
            #moment2_y = np.var(std_velo_dis_y)
            moment3_y = sp.skew(std_velo_dis_y)
            moment4_y = sp.kurtosis(std_velo_dis_y)
            hist_y = plt.hist(std_velo_dis_y, bins=self.normal_bins_size[idx], density=True, log=True)
            plt.close()

            std_velo_dis_z = (velo_z - np.mean(velo_z)) / np.std(velo_z)
            #moment1_z = np.mean(std_velo_dis_z)
            #moment2_z = np.var(std_velo_dis_z)
            moment3_z = sp.skew(std_velo_dis_z)
            moment4_z = sp.kurtosis(std_velo_dis_z)
            hist_z = plt.hist(std_velo_dis_z, bins=self.normal_bins_size[idx], density=True, log=True)
            param = norm.fit(std_velo_dis_z)
            bins_pdf = hist_z[1][:-1] + (0.5 * np.diff(hist_z[1]))
            vals_pdf = norm.pdf(bins_pdf, loc=param[0], scale=param[1] - 0.1)
            #r2_z = r2_score(hist_z[0], vals_pdf)
            plt.close()

            plt.figure(figsize=(4.5,3),dpi=300)
            ax = plt.gca()        
            ax.plot((hist_x[1][:-1] + (0.5 * np.diff(hist_x[1])))[::], (hist_x[0])[::], linestyle='None', markersize=4, marker="o", c='r', mfc="None", alpha=0.5, label="x-direction")
            ax.plot((hist_y[1][:-1] + (0.5 * np.diff(hist_y[1])))[::], (hist_y[0])[::], linestyle='None', markersize=4, marker="^", c='b', mfc="None", alpha=0.5, label="y-direction")
            ax.plot((hist_z[1][:-1] + (0.5 * np.diff(hist_z[1])))[::], (hist_z[0])[::], linestyle='None', markersize=4, marker="s", c='g', mfc="None", alpha=0.5, label="z-direction")
            ax.plot(bins_pdf, vals_pdf, lw=1.0, c='k', ls='-')

            ax.set_yscale('log')
            ax.set_xlim([-5.0,5.0])
            ax.set_xticks([-5.0, -2.5, 0, 2.5, 5.0])
            ax.set_xticklabels([-5.0, -2.5, 0, 2.5, 5.0])
            ax.set_ylim([10**-3,10**-0])
            #ax.set_yticks([0, 0.5, 1.0])
            #ax.set_yticklabels([0, 12.5, 25])
            ax.set_xlabel('normalized velocity')
            ax.set_ylabel('probability density')
            ax.legend(bbox_transform=ax.transAxes, bbox_to_anchor=(0.68, 0.1), loc='lower right', borderaxespad=0)
            ax.text(0.7, 0.7, '$\it{s}_{x}$'+' = {0:.3f} '.format(moment3_x)+ \
                '\n$\it{s}_{y}$'+' = {0:.3f} '.format(moment3_y)+ \
                    '\n$\it{s}_{z}$'+' = {0:.3f} '.format(moment3_z), \
                        rotation=0, transform=ax.transAxes, multialignment='left')
            ax.text(0.05, 0.7, '$\it{k}_{x}$'+' = {0:.3f} '.format(moment4_x)+ \
                '\n$\it{k}_{y}$'+' = {0:.3f} '.format(moment4_y)+ \
                    '\n$\it{k}_{z}$'+' = {0:.3f} '.format(moment4_z), \
                        rotation=0, transform=ax.transAxes, multialignment='left')
            plt.savefig('_post/{0}/{1}_velo_dist_{2:0>5}_ms_.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), velo_measurement), bbox_inches='tight')
            #plt.show()
            plt.close()


    def accumulated_velo_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            accu_velo_x_1 = []
            accu_velo_y_1 = []
            accu_velo_z_1 = []
            accu_velo_x_2 = []
            accu_velo_y_2 = []
            accu_velo_z_2 = []
            count_1 = 0
            count_2 = 0

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                    
                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        accu_velo_x_1.extend(df['velo_x (m/s)'].values)    
                        accu_velo_y_1.extend(df['velo_y (m/s)'].values)
                        accu_velo_z_1.extend(df['velo_z (m/s)'].values)
                        count_1 += 1
                
                    if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                    
                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        accu_velo_x_2.extend(df['velo_x (m/s)'].values)    
                        accu_velo_y_2.extend(df['velo_y (m/s)'].values)
                        accu_velo_z_2.extend(df['velo_z (m/s)'].values)
                        count_2 += 1
                        
            std_velo_dis_x_1 = (np.array(accu_velo_x_1) - np.mean(accu_velo_x_1)) / np.std(accu_velo_x_1)
            #moment1_x = np.mean(std_velo_dis_x_1)
            #moment2_x = np.var(std_velo_dis_x_1)
            moment3_x = sp.skew(std_velo_dis_x_1)
            moment4_x = sp.kurtosis(std_velo_dis_x_1)
            hist_x = plt.hist(std_velo_dis_x_1, bins=self.normal_bins_size[idx], density=True, log=True)
            
            std_velo_dis_y_1 = (np.array(accu_velo_y_1) - np.mean(accu_velo_y_1)) / np.std(accu_velo_y_1)
            #moment1_y = np.mean(std_velo_dis_y_1)
            #moment2_y = np.var(std_velo_dis_y_1)
            moment3_y = sp.skew(std_velo_dis_y_1)
            moment4_y = sp.kurtosis(std_velo_dis_y_1)
            hist_y = plt.hist(std_velo_dis_y_1, bins=self.normal_bins_size[idx], density=True, log=True)
            
            std_velo_dis_z_1 = (np.array(accu_velo_z_1) - np.mean(accu_velo_z_1)) / np.std(accu_velo_z_1)
            #moment1_z = np.mean(std_velo_dis_z_1)
            #moment2_z = np.var(std_velo_dis_z_1)
            moment3_z = sp.skew(std_velo_dis_z_1)
            moment4_z = sp.kurtosis(std_velo_dis_z_1)
            hist_z = plt.hist(std_velo_dis_z_1, bins=self.normal_bins_size[idx], density=True, log=True)
            param = norm.fit(std_velo_dis_z_1)
            bins_pdf = hist_z[1][:-1] + (0.5 * np.diff(hist_z[1]))
            vals_pdf = norm.pdf(bins_pdf, loc=param[0], scale=param[1] - 0.1)
            #r2_z = r2_score(hist_z[0], vals_pdf)
            plt.close()

            plt.figure(figsize=(4.5,3),dpi=300)
            ax = plt.gca()        
            ax.plot((hist_x[1][:-1] + (0.5 * np.diff(hist_x[1])))[::], (hist_x[0])[::], linestyle='None', markersize=4, marker="o", c='r', mfc="None", alpha=0.5, label="x-direction")
            ax.plot((hist_y[1][:-1] + (0.5 * np.diff(hist_y[1])))[::], (hist_y[0])[::], linestyle='None', markersize=4, marker="^", c='b', mfc="None", alpha=0.5, label="y-direction")
            ax.plot((hist_z[1][:-1] + (0.5 * np.diff(hist_z[1])))[::], (hist_z[0])[::], linestyle='None', markersize=4, marker="s", c='g', mfc="None", alpha=0.5, label="z-direction")
            ax.plot(bins_pdf, vals_pdf, lw=1.0, c='k', ls='-')

            ax.set_yscale('log')
            ax.set_xlim([-5.0,5.0])
            ax.set_xticks([-5.0, -2.5, 0, 2.5, 5.0])
            ax.set_xticklabels([-5.0, -2.5, 0, 2.5, 5.0])
            ax.set_ylim([10**-3,10**-0])
            #ax.set_yticks([0, 0.5, 1.0])
            #ax.set_yticklabels([0, 12.5, 25])
            ax.set_xlabel('normalized velocity')
            ax.set_ylabel('probability density')
            ax.legend(bbox_transform=ax.transAxes, bbox_to_anchor=(0.68, 0.1), loc='lower right', borderaxespad=0)
            ax.text(0.7, 0.7, '$\it{s}_{x}$'+' = {0:.3f} '.format(moment3_x)+ \
                '\n$\it{s}_{y}$'+' = {0:.3f} '.format(moment3_y)+ \
                    '\n$\it{s}_{z}$'+' = {0:.3f} '.format(moment3_z), \
                        rotation=0, transform=ax.transAxes, multialignment='left')
            ax.text(0.05, 0.7, '$\it{k}_{x}$'+' = {0:.3f} '.format(moment4_x)+ \
                '\n$\it{k}_{y}$'+' = {0:.3f} '.format(moment4_y)+ \
                    '\n$\it{k}_{z}$'+' = {0:.3f} '.format(moment4_z), \
                        rotation=0, transform=ax.transAxes, multialignment='left')
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_excitation_velo_dist.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_1, count_1 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()


            std_velo_dis_x_2 = (np.array(accu_velo_x_2) - np.mean(accu_velo_x_2)) / np.std(accu_velo_x_2)
            #moment1_x = np.mean(std_velo_dis_x_2)
            #moment2_x = np.var(std_velo_dis_x_2)
            moment3_x = sp.skew(std_velo_dis_x_2)
            moment4_x = sp.kurtosis(std_velo_dis_x_2)
            hist_x = plt.hist(std_velo_dis_x_2, bins=self.normal_bins_size[idx], density=True, log=True)
            
            std_velo_dis_y_2 = (np.array(accu_velo_y_2) - np.mean(accu_velo_y_2)) / np.std(accu_velo_y_2)
            #moment1_y = np.mean(std_velo_dis_y_2)
            #moment2_y = np.var(std_velo_dis_y_2)
            moment3_y = sp.skew(std_velo_dis_y_2)
            moment4_y = sp.kurtosis(std_velo_dis_y_2)
            hist_y = plt.hist(std_velo_dis_y_2, bins=self.normal_bins_size[idx], density=True, log=True)
                    
            std_velo_dis_z_2 = (np.array(accu_velo_z_2) - np.mean(accu_velo_z_2)) / np.std(accu_velo_z_2)
            #moment1_z = np.mean(std_velo_dis_z_2)
            #moment2_z = np.var(std_velo_dis_z_2)
            moment3_z = sp.skew(std_velo_dis_z_2)
            moment4_z = sp.kurtosis(std_velo_dis_z_2)
            hist_z = plt.hist(std_velo_dis_z_2, bins=self.normal_bins_size[idx], density=True, log=True)
            param = norm.fit(std_velo_dis_z_2)
            bins_pdf = hist_z[1][:-1] + (0.5 * np.diff(hist_z[1]))
            vals_pdf = norm.pdf(bins_pdf, loc=param[0], scale=param[1] - 0.03)
            #r2_z = r2_score(hist_z[0], vals_pdf)
            plt.close() 

            plt.figure(figsize=(4.5,3),dpi=300)
            ax = plt.gca()        
            ax.plot((hist_x[1][:-1] + (0.5 * np.diff(hist_x[1])))[::], (hist_x[0])[::], linestyle='None', markersize=4, marker="o", c='r', mfc="None", alpha=0.5, label="x-direction")
            ax.plot((hist_y[1][:-1] + (0.5 * np.diff(hist_y[1])))[::], (hist_y[0])[::], linestyle='None', markersize=4, marker="^", c='b', mfc="None", alpha=0.5, label="y-direction")
            ax.plot((hist_z[1][:-1] + (0.5 * np.diff(hist_z[1])))[::], (hist_z[0])[::], linestyle='None', markersize=4, marker="s", c='g', mfc="None", alpha=0.5, label="z-direction")
            ax.plot(bins_pdf, vals_pdf, lw=1.0, c='k', ls='-')

            ax.set_yscale('log')
            ax.set_xlim([-5.0,5.0])
            ax.set_xticks([-5.0, -2.5, 0, 2.5, 5.0])
            ax.set_xticklabels([-5.0, -2.5, 0, 2.5, 5.0])
            ax.set_ylim([10**-3,10**-0])
            #ax.set_yticks([0, 0.5, 1.0])
            #ax.set_yticklabels([0, 12.5, 25])
            ax.set_xlabel('normalized velocity')
            ax.set_ylabel('probability density')
            ax.legend(bbox_transform=ax.transAxes, bbox_to_anchor=(0.68, 0.1), loc='lower right', borderaxespad=0)
            ax.text(0.7, 0.7, '$\it{s}_{x}$'+' = {0:.3f} '.format(moment3_x)+ \
                '\n$\it{s}_{y}$'+' = {0:.3f} '.format(moment3_y)+ \
                    '\n$\it{s}_{z}$'+' = {0:.3f} '.format(moment3_z), \
                        rotation=0, transform=ax.transAxes, multialignment='left')
            ax.text(0.05, 0.7, '$\it{k}_{x}$'+' = {0:.3f} '.format(moment4_x)+ \
                '\n$\it{k}_{y}$'+' = {0:.3f} '.format(moment4_y)+ \
                    '\n$\it{k}_{z}$'+' = {0:.3f} '.format(moment4_z), \
                        rotation=0, transform=ax.transAxes, multialignment='left')
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_cooling_velo_dist.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_2, count_2 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()


    def velo_dist_3_dir(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            velo_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), velo_measurement))        
            fig, axes = plt.subplots(3, 1, figsize=(4.5,12), dpi=300, sharex=True, sharey=False)
            fig.subplots_adjust(hspace=0.15)

            velo_dis_x = np.where(df['velo_x (m/s)'].values == 0., df['velo_x (m/s)'].values, df['velo_x (m/s)'].values * 1000)            
            hist_x = axes[0].hist(velo_dis_x, bins=self.normal_bins_size[idx], density=True, log=True)
            param = norm.fit(velo_dis_x)
            if ((velo_dis_x.min() != 0.) or (velo_dis_x.max() != 0.)):   
                bins_pdf = hist_x[1][:-1] + (0.5 * np.diff(hist_x[1]))
                vals_pdf = norm.pdf(bins_pdf, loc=param[0], scale=param[1])
                axes[0].plot(bins_pdf, vals_pdf, lw=2.0, c='r', ls='dashed')
                axes[0].text(0.73, 0.9, '$\it{R}$$^2$'+'= {0:.4f}'.format(r2_score(hist_x[0], vals_pdf)), rotation=0, transform=axes[0].transAxes, bbox=dict(lw=0.5, facecolor='w', alpha=1.))
            
            #axes[0].set_xlim(-120,120)
            #axes[0].set_xticks([-120, -80, -40, 0, 40, 80, 120])
            #axes[0].set_xticklabels([-120, -80, -40, 0, 40, 80, 120], rotation=0)
            #axes[0].set_ylim(10**-4,10**-2*3)
            axes[0].set_ylabel('probability density')
            axes[0].set_title('x direction', loc='left')           
            
            velo_dis_y = np.where(df['velo_y (m/s)'].values == 0., df['velo_y (m/s)'].values, df['velo_y (m/s)'].values * 1000)
            hist_y = axes[1].hist(velo_dis_y, bins=self.normal_bins_size[idx], density=True, log=True)
            param = norm.fit(velo_dis_y)
            if ((velo_dis_y.min() != 0.) or (velo_dis_y.max() != 0.)):   
                bins_pdf = hist_y[1][:-1] + (0.5 * np.diff(hist_y[1]))
                vals_pdf = norm.pdf(bins_pdf, loc=param[0], scale=param[1])
                axes[1].plot(bins_pdf, vals_pdf, lw=2.0, c='r', ls='dashed')
                axes[1].text(0.73, 0.9, '$\it{R}$$^2$'+'= {0:.4f}'.format(r2_score(hist_y[0], vals_pdf)), rotation=0, transform=axes[1].transAxes, bbox=dict(lw=0.5, facecolor='w', alpha=1.))
            
            #axes[1].set_xlim(-120,120)
            #axes[1].set_xticks([-120, -80, -40, 0, 40, 80, 120])
            #axes[1].set_xticklabels([-120, -80, -40, 0, 40, 80, 120], rotation=0)
            #axes[1].set_ylim(10**-4,10**-2*3)
            axes[1].set_ylabel('probability density')
            axes[1].set_title('y direction', loc='left')

            velo_dis_z = np.where(df['velo_z (m/s)'].values == 0., df['velo_z (m/s)'].values, df['velo_z (m/s)'].values * 1000)
            hist_z = axes[2].hist(velo_dis_z, bins=self.normal_bins_size[idx], density=True, log=True)
            param = norm.fit(velo_dis_z)
            if ((velo_dis_z.min() != 0.) or (velo_dis_z.max() != 0.)):   
                bins_pdf = hist_z[1][:-1] + (0.5 * np.diff(hist_z[1]))
                vals_pdf = norm.pdf(bins_pdf, loc=param[0], scale=param[1])
                axes[2].plot(bins_pdf, vals_pdf, lw=2.0, c='r', ls='dashed')
                axes[2].text(0.73, 0.9, '$\it{R}$$^2$'+'= {0:.4f}'.format(r2_score(hist_z[0], vals_pdf)), rotation=0, transform=axes[2].transAxes, bbox=dict(lw=0.5, facecolor='w', alpha=1.))
            
            #axes[2].set_xlim(-120,120)
            #axes[2].set_xticks([-120, -80, -40, 0, 40, 80, 120])
            #axes[2].set_xticklabels([-120, -80, -40, 0, 40, 80, 120], rotation=0)
            #axes[2].set_ylim(10**-4,10**-2*3)
            axes[2].set_xlabel('velocity  $\it{mm/s}$')
            axes[2].set_ylabel('probability density')
            axes[2].set_title('z direction', loc='left')

            plt.savefig('_post/{0}/{1}_velo_dist_3_dir_{2:0>5}_ms_.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), velo_measurement), bbox_inches='tight')
            plt.close()


    def maxwellian_velo_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            velo_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), velo_measurement))        
            std_velo_dis = df['velo_ab (m/s)'].values / np.mean(df['velo_ab (m/s)'].values)
            hist = plt.hist(std_velo_dis, bins=self.normal_bins_size[idx], density=True, log=True)
            bins_pdf = hist[1][:-1] + (0.5 * np.diff(hist[1]))
            param = maxwell.fit(std_velo_dis)
            vals_pdf = maxwell.pdf(bins_pdf, loc=param[0], scale=param[1])
            #r2_2 = r2_score(hist_2[0], vals_pdf)
            exp_x = np.arange(2.5, 5.1, 0.1)
            num_A = 0.52
            num_B = 2.16
            plt.close()  
            
            plt.figure(figsize=(4.5,3),dpi=300)
            ax = plt.gca()        
            ax.plot((hist[1][:-1] + (0.5 * np.diff(hist[1])))[::], (hist[0])[::], linestyle='None', markersize=4, marker="o", c='r', mfc="None", alpha=0.5)
            ax.plot(bins_pdf, vals_pdf, lw=1.0, c='K', ls='-')  
            ax.plot(exp_x, np.exp(-num_A*(exp_x**num_B)), lw=2.0, c='r', ls='-') 
            
            ax.set_yscale('log')
            ax.set_xlim([0.0,5.0])
            ax.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
            ax.set_xticklabels([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
            ax.set_ylim([10**-3,10**-0])
            #ax.set_yticks([0, 0.5, 1.0])
            #ax.set_yticklabels([0, 12.5, 25])
            ax.set_xlabel('velocity / average velocity $\it{c}$')
            ax.set_ylabel('probability density')
            ax.text(0.28, 0.15, 'Maxwell- \nBoltzmann', rotation=0, multialignment='left', \
                verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes)
            arrow_dict = dict(color='r', width=0.1, headwidth=5, headlength=5)
            ax.annotate('~ $\it{e}^{\it{-\u03b1c^{\it{\u03b2}}}}}$ \n'+'$(\it{\u03b2}$'+' = {0:.2f}) '.format(num_B), \
                            xy = (2.9, 2*10**-2), xytext = (3.8, 6*10**-2), color = 'r', arrowprops = arrow_dict)
            plt.savefig('_post/{0}/{1}_velo_dist_maxwell_{2:0>5}_ms_.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), velo_measurement), bbox_inches='tight')        
            #plt.show()
            plt.close()


    def accumulated_maxwellian_velo_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):
    
            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            accu_velo_1 = []
            accu_velo_2 = []
            count_1 = 0
            count_2 = 0

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                    
                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        accu_velo_1.extend(df['velo_ab (m/s)'].values)    
                        count_1 += 1
                
                    if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                    
                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        accu_velo_2.extend(df['velo_ab (m/s)'].values)    
                        count_2 += 1
    
            std_velo_dis_1 = np.array(accu_velo_1) / np.mean(accu_velo_1)
            hist_1 = plt.hist(std_velo_dis_1, bins=self.normal_bins_size[idx], density=True, log=True)
            num_A_1 = 1.34
            num_B_1 = 1.05
            plt.close()

            std_velo_dis_2 = np.array(accu_velo_2) / np.mean(accu_velo_2)
            hist_2 = plt.hist(std_velo_dis_2, bins=self.normal_bins_size[idx], density=True, log=True)
            bins_pdf = hist_2[1][:-1] + (0.5 * np.diff(hist_2[1]))
            param = maxwell.fit(std_velo_dis_2)
            vals_pdf = maxwell.pdf(bins_pdf, loc=param[0], scale=param[1])
            #r2_2 = r2_score(hist_2[0], vals_pdf)
            exp_x = np.arange(2.5, 5.1, 0.1)
            num_A_2 = 0.52
            num_B_2 = 2.16
            plt.close()  
            

            plt.figure(figsize=(4.5,3),dpi=300)
            ax = plt.gca()        
            ax.plot((hist_1[1][:-1] + (0.5 * np.diff(hist_1[1])))[::], (hist_1[0])[::], linestyle='None', markersize=4, marker="o", c='r', mfc="None", alpha=0.5, label="end of turning-on")
            ax.plot((hist_2[1][:-1] + (0.5 * np.diff(hist_2[1])))[::], (hist_2[0])[::], linestyle='None', markersize=4, marker="^", c='b', mfc="None", alpha=0.5, label="end of turning-off")
            ax.plot(bins_pdf, vals_pdf, lw=1.0, c='K', ls='-')  
            ax.plot(exp_x, np.exp(-num_A_1*(exp_x**num_B_1)), lw=2.0, c='r', ls='-') 
            ax.plot(exp_x, np.exp(-num_A_2*(exp_x**num_B_2)), lw=2.0, c='b', ls='-')

            ax.set_yscale('log')
            ax.set_xlim([0.0,5.0])
            ax.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
            ax.set_xticklabels([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
            ax.set_ylim([10**-3,10**-0])
            #ax.set_yticks([0, 0.5, 1.0])
            #ax.set_yticklabels([0, 12.5, 25])
            ax.set_xlabel('velocity / average velocity $\it{c}$')
            ax.set_ylabel('probability density')
            ax.legend(bbox_transform=ax.transAxes, bbox_to_anchor=(0.99, 0.8), loc='lower right', borderaxespad=0)
            ax.text(0.28, 0.15, 'Maxwell- \nBoltzmann', rotation=0, multialignment='left', \
                verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes)

            arrow_dict_1 = dict(color='r', width=0.1, headwidth=5, headlength=5)
            ax.annotate('~ $\it{e}^{\it{-\u03b1c^{\it{\u03b2}}}}}$ \n'+'$(\it{\u03b2}$'+' = {0:.2f}) '.format(num_B_1), \
                            xy = (2.9, 2*10**-2), xytext = (3.8, 6*10**-2), color = 'r', arrowprops = arrow_dict_1)
            
            arrow_dict_2 = dict(color='b', width=0.1, headwidth=5, headlength=5)
            ax.annotate('~ $\it{e}^{\it{-\u03b1c^{\it{\u03b2}}}}}$ \n'+'$(\it{\u03b2}$'+' = {0:.2f}) '.format(num_B_2), \
                            xy = (3.2, 2*10**-3), xytext = (3.8, 1.2*10**-2), color = 'b', arrowprops = arrow_dict_2)

            #ax.text(0.7, 0.4, '~ $\it{e}^{\it{-\u03b1c^{\it{\u03b2}}}}}$ \n'+'$(\it{\u03b2}$'+' = {0:.2f}) '.format(num_B_1), \
            #            rotation=0, verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes, multialignment='left')
            #ax.text(0.7, 0.7, '~ $\it{e}^{\it{-\u03b1c^{\it{\u03b2}}}}}$ \n'+'$(\it{\u03b2}$'+' = {0:.2f}) '.format(num_B_2), \
            #            rotation=0, verticalalignment='bottom', horizontalalignment='left', transform=ax.transAxes, multialignment='left')
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_cooling_velo_dist_maxwell.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_1, count_1 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()


    def velo_dir(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):
            
            velo_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), velo_measurement))   

            velo_x_unit = np.divide(df['velo_x (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)
            velo_y_unit = np.divide(df['velo_y (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)
            velo_z_unit = np.divide(df['velo_z (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)
            
            x_dir = np.array(plt.hist(velo_x_unit, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            y_dir = np.array(plt.hist(velo_y_unit, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            z_dir = np.array(plt.hist(velo_z_unit, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            plt.close()

            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(x_dir[1][:-1] + (0.5 * np.diff(x_dir[1])), x_dir[0], c=temp_cmap(1.0), ls='-', lw = 1.25, label='velocity x-direction')
            plt.plot(y_dir[1][:-1] + (0.5 * np.diff(y_dir[1])), y_dir[0], c=temp_cmap(0.8), ls='--', lw = 2.0, label='velocity y-direction')
            plt.plot(z_dir[1][:-1] + (0.5 * np.diff(z_dir[1])), z_dir[0], c=temp_cmap(0.6), ls=':', lw = 1.0, label='velocity z-direction')
            plt.xlabel('direction')
            plt.ylabel('probability density')
            plt.xlim(xmin=-1.0,xmax=1.0)
            plt.xticks([-1.0, -0.5, 0, 0.5, 1.0], [-1.0, -0.5, 0, 0.5, 1.0],rotation=0)
            plt.ylim(ymin=0,ymax=1.0)
            plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 0.2, 0.4, 0.6, 0.8, 1.0],rotation=0)
            plt.legend(bbox_to_anchor=(0.98, 0.72), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_velo_dir_{2}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), velo_measurement),bbox_inches='tight')
            plt.close()


    def accumulated_velo_dir(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count_1 = 0
            count_2 = 0
            accu_velo_x_unit_1 = []
            accu_velo_y_unit_1 = []
            accu_velo_z_unit_1 = []
            accu_velo_x_unit_2 = []
            accu_velo_y_unit_2 = []
            accu_velo_z_unit_2 = []
            
            
            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    
                    if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                    
                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        accu_velo_x_unit_1.extend((np.divide(df['velo_x (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)).tolist())    
                        accu_velo_y_unit_1.extend((np.divide(df['velo_y (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)).tolist())
                        accu_velo_z_unit_1.extend((np.divide(df['velo_z (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)).tolist())
                        count_1 += 1
                
                    if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                    
                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))
                        accu_velo_x_unit_2.extend((np.divide(df['velo_x (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)).tolist())    
                        accu_velo_y_unit_2.extend((np.divide(df['velo_y (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)).tolist())
                        accu_velo_z_unit_2.extend((np.divide(df['velo_z (m/s)'].values, df['velo_ab (m/s)'].values, out=np.zeros(len(df)), where=df['velo_ab (m/s)'].values!=0.)).tolist())
                        count_2 += 1

            x_dir_1 = np.array(plt.hist(accu_velo_x_unit_1, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            y_dir_1 = np.array(plt.hist(accu_velo_y_unit_1, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            z_dir_1 = np.array(plt.hist(accu_velo_z_unit_1, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            x_dir_2 = np.array(plt.hist(accu_velo_x_unit_2, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            y_dir_2 = np.array(plt.hist(accu_velo_y_unit_2, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            z_dir_2 = np.array(plt.hist(accu_velo_z_unit_2, bins=int(self.normal_bins_size[idx]*0.5), density=True))
            plt.close()


            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(x_dir_1[1][:-1] + (0.5 * np.diff(x_dir_1[1])), x_dir_1[0], c=temp_cmap(1.0), ls='-', lw = 1.25, label='velocity x-direction')
            plt.plot(y_dir_1[1][:-1] + (0.5 * np.diff(y_dir_1[1])), y_dir_1[0], c=temp_cmap(0.8), ls='--', lw = 2.0, label='velocity y-direction')
            plt.plot(z_dir_1[1][:-1] + (0.5 * np.diff(z_dir_1[1])), z_dir_1[0], c=temp_cmap(0.6), ls=':', lw = 1.0, label='velocity z-direction')
            plt.xlabel('direction')
            plt.ylabel('probability density')
            plt.xlim(xmin=-1.0,xmax=1.0)
            plt.xticks([-1.0, -0.5, 0, 0.5, 1.0], [-1.0, -0.5, 0, 0.5, 1.0],rotation=0)
            plt.ylim(ymin=0,ymax=1.0)
            plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 0.2, 0.4, 0.6, 0.8, 1.0],rotation=0)
            plt.legend(bbox_to_anchor=(0.98, 0.72), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_excitation_velo_dir.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_1, count_1 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()
                        
            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(x_dir_2[1][:-1] + (0.5 * np.diff(x_dir_2[1])), x_dir_2[0], c=temp_cmap(1.0), ls='-', lw = 1.25, label='velocity x-direction')
            plt.plot(y_dir_2[1][:-1] + (0.5 * np.diff(y_dir_2[1])), y_dir_2[0], c=temp_cmap(0.8), ls='--', lw = 2.0, label='velocity y-direction')
            plt.plot(z_dir_2[1][:-1] + (0.5 * np.diff(z_dir_2[1])), z_dir_2[0], c=temp_cmap(0.6), ls=':', lw = 1.0, label='velocity z-direction')
            plt.xlabel('direction')
            plt.ylabel('probability density')
            plt.xlim(xmin=-1.0,xmax=1.0)
            plt.xticks([-1.0, -0.5, 0, 0.5, 1.0], [-1.0, -0.5, 0, 0.5, 1.0],rotation=0)
            plt.ylim(ymin=0,ymax=1.0)
            plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], [0, 0.2, 0.4, 0.6, 0.8, 1.0],rotation=0)
            plt.legend(bbox_to_anchor=(0.98, 0.72), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_cooling_velo_dir.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_2, count_2 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()
            

    def cell_layer_velo_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):
            
            velo_dist_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), velo_dist_measurement))   
            df_pos = pd.DataFrame()
            
            if(self.cubic_cell_bool[idx]):
                half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
            
            else:
                sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03   # adjust for calculations condition
                

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

            else:
                rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)
                ave_velo_in_each_cell = np.zeros(sphere_radius)
                ptcl_num_in_each_cell = np.zeros(sphere_radius)
                

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
            plt.figure(figsize=(4.5,3),dpi=300)    

            if(self.cubic_cell_bool[idx]):
                plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell, c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('cell layer')
            else:
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell, c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('layer')            
            plt.ylabel('velocity  $\it{mm/s}$')        
            plt.xlim(xmin=5,xmax=25)
            plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
            #plt.ylim(ymin=0,ymax=45.0)
            #plt.yticks([0, 15, 30, 45], [0, 15, 30, 45],rotation=0)
            
            plt.savefig('_post/{0}/{1}_velo_cell_{2}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), velo_dist_measurement), bbox_inches='tight')
            #plt.show()
            plt.close()


    def accumulated_cell_layer_velo_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count_1 = 0
            count_2 = 0

            if(self.cubic_cell_bool[idx]):

                half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
                ave_velo_in_each_cell_1 = np.zeros(half_cubic_length)
                ptcl_num_in_each_cell_1 = np.zeros(half_cubic_length)
                ave_velo_in_each_cell_2 = np.zeros(half_cubic_length)
                ptcl_num_in_each_cell_2 = np.zeros(half_cubic_length)

            else:
                sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
                ave_velo_in_each_cell_1 = np.zeros(sphere_radius)
                ptcl_num_in_each_cell_1 = np.zeros(sphere_radius)
                ave_velo_in_each_cell_2 = np.zeros(sphere_radius)
                ptcl_num_in_each_cell_2 = np.zeros(sphere_radius)
        

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                file_time = re.findall(r'\d+', str(val_csv.name))[0]

                if (idx_csv != 0) and (int(self.time_start_accum[idx] * 1000) < int(file_time)) and (int(self.end_excitation[idx] * 1000) >= int(file_time)):
                    if (((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0) or ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0)):
                    
                        df = pd.read_csv('{0}/OUTPUT_2/{1}'.format(str(self.tgt_dirs[idx]), val_csv.name))                        
                        df_pos = pd.DataFrame()
            
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


                        if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                            count_1 += 1
                        if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                            count_2 += 1  
                            
                        for idx_ptcl in range (len(df_pos)):

                            if(self.cubic_cell_bool[idx]):

                                temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                                dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<=half_cubic_length - 1) else (temp_pos[0] - half_cubic_length).astype(int)
                                dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<=half_cubic_length - 1) else (temp_pos[1] - half_cubic_length).astype(int)
                                dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<=half_cubic_length - 1) else (temp_pos[2] - half_cubic_length).astype(int)
                                in_cell = max([dis_x, dis_y, dis_z])

                                if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                        
                                    ave_velo_in_each_cell_1[in_cell] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                    ptcl_num_in_each_cell_1[in_cell] += 1  
                                    
                                if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                
                                    ave_velo_in_each_cell_2[in_cell] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                    ptcl_num_in_each_cell_2[in_cell] += 1  
                                    
                            else:
                                if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                        
                                    ave_velo_in_each_cell_1[rel_pos[idx_ptcl]] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                    ptcl_num_in_each_cell_1[rel_pos[idx_ptcl]] += 1
                                    
                                if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                
                                    ave_velo_in_each_cell_2[rel_pos[idx_ptcl]] += df.at[idx_ptcl, 'velo_ab (m/s)']
                                    ptcl_num_in_each_cell_2[rel_pos[idx_ptcl]] += 1
                                    

            for idx_cell in range (4):
                ave_velo_in_each_cell_1[4] += ave_velo_in_each_cell_1[idx_cell]
                ptcl_num_in_each_cell_1[4] += ptcl_num_in_each_cell_1[idx_cell]
                ave_velo_in_each_cell_1[idx_cell] = 0
                ptcl_num_in_each_cell_1[idx_cell] = 0
                ave_velo_in_each_cell_2[4] += ave_velo_in_each_cell_2[idx_cell]
                ptcl_num_in_each_cell_2[4] += ptcl_num_in_each_cell_2[idx_cell]
                ave_velo_in_each_cell_2[idx_cell] = 0
                ptcl_num_in_each_cell_2[idx_cell] = 0
            
            output_ave_velo_in_each_cell_1 = 1000 * np.divide(ave_velo_in_each_cell_1, ptcl_num_in_each_cell_1, out=np.zeros_like(ave_velo_in_each_cell_1), where=ptcl_num_in_each_cell_1!=0)
            output_ave_velo_in_each_cell_2 = 1000 * np.divide(ave_velo_in_each_cell_2, ptcl_num_in_each_cell_2, out=np.zeros_like(ave_velo_in_each_cell_2), where=ptcl_num_in_each_cell_2!=0)
            
            plt.figure(figsize=(4.5,3),dpi=300)
            if(self.cubic_cell_bool[idx]):
                plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell_1, c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('cell layer')
            else:
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell_1, c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('layer')        
            plt.ylabel('velocity  $\it{mm/s}$')       
            plt.xlim(xmin=5,xmax=25)
            plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
            #plt.ylim(ymin=0,ymax=50.0)
            #plt.yticks([0, 15, 30, 45], [0, 15, 30, 45],rotation=0)
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_excitation_velo_cell.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_1, count_1 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()


            plt.figure(figsize=(4.5,3),dpi=300)
            if(self.cubic_cell_bool[idx]):
                plt.plot(np.array(range(1,half_cubic_length+1)), output_ave_velo_in_each_cell_2, c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('cell layer')
            else:
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_velo_in_each_cell_2, c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('layer')        
            plt.ylabel('velocity  $\it{mm/s}$')      
            plt.xlim(xmin=5,xmax=25)
            plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
            #plt.ylim(ymin=0,ymax=50.0)
            #plt.yticks([0, 15, 30, 45], [0, 15, 30, 45],rotation=0)
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_cooling_velo_cell.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_2, count_2 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()
            

    

            