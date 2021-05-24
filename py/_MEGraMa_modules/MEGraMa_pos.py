import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib as path
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

temp_cmap = plt.get_cmap("Blues")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams["mathtext.fontset"] = 'dejavusans'
plt.rcParams['font.size'] = 10
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
cr_path = path.Path()


class Pos_data_output:

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

    
    def pos_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):
    
            # read csv file
            pos_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), pos_measurement)) 
            df_pos = pd.DataFrame()
            cubic_length = int((self.cell_length[idx] * 1000) + 0.5)
            sphere_cell_radius = self.cell_length[idx] * 0.5
            division_number_to_radius = 30
            division_number_to_radian = 60
            divided_r = np.array([ np.cbrt(idx_r+1) for idx_r in range(division_number_to_radius)]) / np.max(np.array([ np.cbrt(idx_r+1) for idx_r in range(division_number_to_radius)]))
            vals_r_base = np.append(0, divided_r[:-1])
            vals_r_norm = divided_r - vals_r_base
            vals_r_square = ((vals_r_base+vals_r_norm) * sphere_cell_radius) * ((vals_r_base+vals_r_norm) * sphere_cell_radius)
            vals_rad_base = np.arccos((1. - np.linspace(0, 2, division_number_to_radian+1)[:-1]))
            vals_rad_norm = np.append(vals_rad_base[1:], np.pi) - vals_rad_base
            vals_rad_sum = vals_rad_base + vals_rad_norm
            
            if(self.cubic_cell_bool[idx]):
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
            else:
                df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03

            if((self.vibration_on[idx]) and (pos_measurement <= int((self.end_excitation[idx] * 1000) + 0.5))):      

                temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                vib_amplitude_x = temp_amplitude[idx] * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                vib_amplitude_y = temp_amplitude[idx] * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(pos_measurement) * 0.001)
                df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(pos_measurement) * 0.001)
                df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(pos_measurement) * 0.001)

            
            if(self.cubic_cell_bool[idx]):

                df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                ptcl_num_cubic = np.zeros((cubic_length, cubic_length)).astype(int)

                # particle spatial distribution at each time step
                for idx_ptcl in range (len(df_pos)):
                    ptcl_num_cubic[df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']] += 1

                plt.figure(figsize=(2.5,2.0),dpi=300)
                ax = sns.heatmap((ptcl_num_cubic / float(np.max(ptcl_num_cubic))).T, cmap='Blues', vmin=0.0, vmax=1.0) 
                plt.xlabel('y  $\it{mm}$')
                plt.ylabel('z  $\it{mm}$')
                plt.xlim([0,50])
                plt.xticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)
                plt.ylim([0,50])
                plt.yticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)
                
                cbar = ax.collections[0].colorbar
                cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                cbar.set_ticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                cbar.update_ticks() 
                plt.savefig('_post/{0}/{1}_pos_prob_{2}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), pos_measurement),bbox_inches='tight')
                #plt.show()
                plt.close()

            else:
                rel_pos_square = (df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)'])
                rel_pos_unit_x = df_pos['pos_x (m)'] / np.sqrt(rel_pos_square)
                ptcl_num_sphere = np.zeros((division_number_to_radius, division_number_to_radian)).astype(int)

                # particle spatial distribution at each time step
                for idx_ptcl in range (len(df_pos)):

                    r_num = np.where(rel_pos_square[idx_ptcl] <= vals_r_square)[0][0]
                    rad_num = np.where(np.arccos(np.clip(((rel_pos_unit_x[idx_ptcl] * 1.0) + (0 * 0) + (0 * 0)), -1, 1)) <= vals_rad_sum)[0][0]
                    ptcl_num_sphere[r_num][rad_num] += 1.0

                fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(3.0,3.0), dpi=300)
                normed_ptcl_num_sphere = ptcl_num_sphere / float(np.max(ptcl_num_sphere))
                colors = temp_cmap(normed_ptcl_num_sphere)

                for idx_r in range(division_number_to_radius):

                    ax.bar(x=vals_rad_base,
                        width=vals_rad_norm, bottom=vals_r_base[idx_r], height=vals_r_norm[idx_r],
                        color=colors[idx_r], edgecolor='w', linewidth=0.0, align="edge")

                ax.set_xlim([0, np.pi])
                ax.set_xticks([0, np.pi*1/4, np.pi*2/4, np.pi*3/4, np.pi*4/4])
                ax.set_xticklabels([0, r'$\frac{1}{4}$'+'$\it{\u03c0}$', r'$\frac{1}{2}$'+'$\it{\u03c0}$', r'$\frac{3}{4}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
                ax.text(-0.1, 0.4, 'Radian  $\it{\u03b8}$', rotation=90, transform=ax.transAxes)
                ax.set_ylim([0, 1.0])
                ax.set_yticks([0, 0.5, 1.0])
                ax.set_yticklabels([0, 12.5, 25])
                ax.text(0.5, 0.05, 'Radius  $\it{mm}$', rotation=0, transform=ax.transAxes)
                ax.grid(False)
                ax.tick_params(which='both', left='on', top='on', direction='in')

                normed_data = Normalize(vmin=np.min(normed_ptcl_num_sphere), vmax=np.max(normed_ptcl_num_sphere))
                mappable = ScalarMappable(cmap=temp_cmap, norm=normed_data)
                mappable._A = []
                cbar = plt.colorbar(mappable, ax=ax, pad=0.12, shrink=0.5)
                cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                cbar.ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                cbar.outline.set_visible(False)

                plt.box(on=None)
                plt.savefig('_post/{0}/{1}_pos_prob_{2}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), pos_measurement),bbox_inches='tight')
                #plt.show()
                plt.close()


    def accumulated_pos_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            cubic_length = int((self.cell_length[idx] * 1000) + 0.5)            
            sphere_cell_radius = self.cell_length[idx] * 0.5
            division_number_to_radius = 30
            division_number_to_radian = 60
            divided_r = np.array([ np.cbrt(idx_r+1) for idx_r in range(division_number_to_radius)]) / np.max(np.array([ np.cbrt(idx_r+1) for idx_r in range(division_number_to_radius)]))
            vals_r_base = np.append(0, divided_r[:-1])
            vals_r_norm = divided_r - vals_r_base
            vals_r_square = ((vals_r_base+vals_r_norm) * sphere_cell_radius) * ((vals_r_base+vals_r_norm) * sphere_cell_radius)
            vals_rad_base = np.arccos((1. - np.linspace(0, 2, division_number_to_radian+1)[:-1]))
            vals_rad_norm = np.append(vals_rad_base[1:], np.pi) - vals_rad_base
            vals_rad_sum = vals_rad_base + vals_rad_norm
            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count_1 = 0
            count_2 = 0
    
            if(self.cubic_cell_bool[idx]):
                ptcl_num_cubic_1 = np.zeros((cubic_length, cubic_length)).astype(int)
                ptcl_num_cubic_2 = np.zeros((cubic_length, cubic_length)).astype(int)
            else:
                ptcl_num_sphere_1 = np.zeros((division_number_to_radius, division_number_to_radian)).astype(int)
                ptcl_num_sphere_2 = np.zeros((division_number_to_radius, division_number_to_radian)).astype(int)

                
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
                            vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy * (np.pi / 180.))
                            vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy * (np.pi / 180.))
                            vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))    
                            df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                            df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(file_time) * 0.001)
                
                        if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                            count_1 += 1
                        if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                            count_2 += 1  

                        if(self.cubic_cell_bool[idx]):
                            
                            df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                            ptcl_num_cubic_temp = np.zeros((cubic_length, cubic_length)).astype(int)
                            
                            for idx_ptcl in range (len(df_pos)):
                                ptcl_num_cubic_temp[df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']] += 1
        
                            if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                ptcl_num_cubic_1 += ptcl_num_cubic_temp
                            if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                ptcl_num_cubic_2 += ptcl_num_cubic_temp 
        
                        else:
                            rel_pos_square = (df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)'])
                            rel_pos_unit_x = df_pos['pos_x (m)'] / np.sqrt(rel_pos_square)
                            ptcl_num_sphere_temp = np.zeros((division_number_to_radius, division_number_to_radian)).astype(int)
        
                            # particle spatial distribution at each time step
                            for idx_ptcl in range (len(df_pos)):
        
                                r_num = np.where(rel_pos_square[idx_ptcl] <= vals_r_square)[0][0]
                                rad_num = np.where(np.arccos(np.clip(((rel_pos_unit_x[idx_ptcl] * 1.0) + (0 * 0) + (0 * 0)), -1, 1)) <= vals_rad_sum)[0][0]
                                ptcl_num_sphere_temp[r_num][rad_num] += 1
                    
                            if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                ptcl_num_sphere_1 += ptcl_num_sphere_temp
                            if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                ptcl_num_sphere_2 += ptcl_num_sphere_temp

            if(self.cubic_cell_bool[idx]): 

                plt.figure(figsize=(2.5,2.0),dpi=300)
                ax = sns.heatmap((ptcl_num_cubic_1 / float(np.max(ptcl_num_cubic_1))).T, cmap='Blues', vmin=0.0, vmax=1.0) 
                plt.xlabel('y  $\it{mm}$')
                plt.ylabel('z  $\it{mm}$')
                plt.xlim([0,50])
                plt.xticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)
                plt.ylim([0,50])
                plt.yticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)            
                cbar = ax.collections[0].colorbar
                cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                cbar.set_ticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                cbar.update_ticks()
                plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_excitation_pos_prob.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_1, count_1 * len(df)), bbox_inches='tight')
                #plt.show()
                plt.close()

                plt.figure(figsize=(2.5,2.0),dpi=300)
                ax = sns.heatmap((ptcl_num_cubic_2 / float(np.max(ptcl_num_cubic_2))).T, cmap='Blues', vmin=0.0, vmax=1.0) 
                plt.xlabel('y  $\it{mm}$')
                plt.ylabel('z  $\it{mm}$')
                plt.xlim([0,50])
                plt.xticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)
                plt.ylim([0,50])
                plt.yticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)            
                cbar = ax.collections[0].colorbar
                cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                cbar.set_ticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                cbar.update_ticks()
                plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_cooling_pos_prob.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_2, count_2 * len(df)), bbox_inches='tight')
                #plt.show()
                plt.close()

            else:
                        
                fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(3.0,3.0), dpi=300)
                normed_ptcl_num_sphere_1 = ptcl_num_sphere_1 / float(np.max(ptcl_num_sphere_1))
                ave_colors_1 = temp_cmap(normed_ptcl_num_sphere_1)

                for idx_r in range(division_number_to_radius):

                    ax.bar(x=vals_rad_base,
                        width=vals_rad_norm, bottom=vals_r_base[idx_r], height=vals_r_norm[idx_r],
                        color=ave_colors_1[idx_r], edgecolor='w', linewidth=0.0, align="edge")

                ax.set_xlim([0, np.pi])
                ax.set_xticks([0, np.pi*1/4, np.pi*2/4, np.pi*3/4, np.pi*4/4])
                ax.set_xticklabels([0, r'$\frac{1}{4}$'+'$\it{\u03c0}$', r'$\frac{1}{2}$'+'$\it{\u03c0}$', r'$\frac{3}{4}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
                ax.text(-0.1, 0.4, 'Radian  $\it{\u03b8}$', rotation=90, transform=ax.transAxes)
                ax.set_ylim([0, 1.0])
                ax.set_yticks([0, 0.5, 1.0])
                ax.set_yticklabels([0, 12.5, 25])
                ax.text(0.5, 0.05, 'Radius  $\it{mm}$', rotation=0, transform=ax.transAxes)
                ax.grid(False)
                ax.tick_params(which='both', left='on', top='on', direction='in')

                normed_data_1 = Normalize(vmin=np.min(normed_ptcl_num_sphere_1), vmax=np.max(normed_ptcl_num_sphere_1))
                mappable_1 = ScalarMappable(cmap=temp_cmap, norm=normed_data_1)
                mappable_1._A = []
                cbar = plt.colorbar(mappable_1, ax=ax, pad=0.12, shrink=0.5)
                cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                cbar.ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                cbar.outline.set_visible(False)
                plt.box(on=None)
                plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_excitation_pos_prob.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_1, count_1 * len(df)), bbox_inches='tight')
                #plt.show()
                plt.close()


                fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(3.0,3.0), dpi=300)
                normed_ptcl_num_sphere_2 = ptcl_num_sphere_2 / float(np.max(ptcl_num_sphere_2))
                ave_colors_2 = temp_cmap(normed_ptcl_num_sphere_2)

                for idx_r in range(division_number_to_radius):

                    ax.bar(x=vals_rad_base,
                        width=vals_rad_norm, bottom=vals_r_base[idx_r], height=vals_r_norm[idx_r],
                        color=ave_colors_2[idx_r], edgecolor='w', linewidth=0.0, align="edge")

                ax.set_xlim([0, np.pi])
                ax.set_xticks([0, np.pi*1/4, np.pi*2/4, np.pi*3/4, np.pi*4/4])
                ax.set_xticklabels([0, r'$\frac{1}{4}$'+'$\it{\u03c0}$', r'$\frac{1}{2}$'+'$\it{\u03c0}$', r'$\frac{3}{4}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
                ax.text(-0.1, 0.4, 'Radian  $\it{\u03b8}$', rotation=90, transform=ax.transAxes)
                ax.set_ylim([0, 1.0])
                ax.set_yticks([0, 0.5, 1.0])
                ax.set_yticklabels([0, 12.5, 25])
                ax.text(0.5, 0.05, 'Radius  $\it{mm}$', rotation=0, transform=ax.transAxes)
                ax.grid(False)
                ax.tick_params(which='both', left='on', top='on', direction='in')

                normed_data_2 = Normalize(vmin=np.min(normed_ptcl_num_sphere_2), vmax=np.max(normed_ptcl_num_sphere_2))
                mappable_2 = ScalarMappable(cmap=temp_cmap, norm=normed_data_2)
                mappable_2._A = []
                cbar = plt.colorbar(mappable_2, ax=ax, pad=0.12, shrink=0.5)
                cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                cbar.ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                cbar.outline.set_visible(False)
                plt.box(on=None)
                plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_cooling_pos_prob.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_2, count_2 * len(df)), bbox_inches='tight')
                #plt.show()
                plt.close()


    def pos_dist_movie(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            cubic_length = int((self.cell_length[idx] * 1000) + 0.5)
            sphere_cell_radius = self.cell_length[idx] * 0.5
            division_number_to_radius = 30
            division_number_to_radian = 60
            divided_r = np.array([ np.cbrt(idx_r+1) for idx_r in range(division_number_to_radius)]) / np.max(np.array([ np.cbrt(idx_r+1) for idx_r in range(division_number_to_radius)]))
            vals_r_base = np.append(0, divided_r[:-1])
            vals_r_norm = divided_r - vals_r_base
            vals_r_square = ((vals_r_base+vals_r_norm) * sphere_cell_radius) * ((vals_r_base+vals_r_norm) * sphere_cell_radius)
            vals_rad_base = np.arccos((1. - np.linspace(0, 2, division_number_to_radian+1)[:-1]))
            vals_rad_norm = np.append(vals_rad_base[1:], np.pi) - vals_rad_base
            vals_rad_sum = vals_rad_base + vals_rad_norm
            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            dir_path = path.Path(cr_path / '_post' / self.tgt_dirs[idx] / 'pos_pic_file')
            dir_path.mkdir(exist_ok=True)

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                # read csv file
                file_time = re.findall(r'\d+', str(val_csv.name))[0]
                df = pd.read_csv(val_csv)
                df_pos = pd.DataFrame()

                if (self.cubic_cell_bool[idx]):
                    df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.005  # adjust for calculations condition
                else:
                    df_pos[['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] = df.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] - 0.03
        
                if((self.vibration_on[idx]) and ((int(file_time)) <= int((self.end_excitation[idx] * 1000) + 0.5))):

                    temp_amplitude = self.vib_amplitude[idx] * np.cos(self.angle_z[idx] * (np.pi / 180.))
                    vib_amplitude_x = temp_amplitude * np.cos(self.angle_xy[idx] * (np.pi / 180.))
                    vib_amplitude_y = temp_amplitude * np.sin(self.angle_xy[idx] * (np.pi / 180.))
                    vib_amplitude_z = self.vib_amplitude[idx] * np.sin(self.angle_z[idx] * (np.pi / 180.))
                    df_pos['pos_x (m)'] -= vib_amplitude_x * np.sin(2. * np.pi * self.vib_frequency[idx] * float(int(file_time)) * 0.001)
                    df_pos['pos_y (m)'] -= vib_amplitude_y * np.sin(2. * np.pi * self.vib_frequency[idx] * float(int(file_time)) * 0.001)
                    df_pos['pos_z (m)'] -= vib_amplitude_z * np.sin(2. * np.pi * self.vib_frequency[idx] * float(int(file_time)) * 0.001)

                if (self.cubic_cell_bool[idx]):

                    df_pos[['pos_x (mm)', 'pos_y (mm)', 'pos_z (mm)']] = (df_pos.loc[:,['pos_x (m)', 'pos_y (m)', 'pos_z (m)']] * 1000).astype(int)  # from m to mm
                    ptcl_pos_array_2d = np.zeros((cubic_length, cubic_length)).astype(int)

                    # particle spatial distribution at each time step
                    for idx_ptcl in range (len(df_pos)):
                        ptcl_pos_array_2d[df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']] += 1

                    plt.figure(figsize=(2.5,2.0),dpi=300)
                    ax = sns.heatmap((ptcl_pos_array_2d / float(np.max(ptcl_pos_array_2d))).T, cmap='Blues', vmin=0.0, vmax=1.0) 
                    plt.xlabel('y  $\it{mm}$')
                    plt.ylabel('z  $\it{mm}$')
                    plt.xlim([0,50])
                    plt.xticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)
                    plt.ylim([0,50])
                    plt.yticks([0,10,20,30,40,50], [0,10,20,30,40,50],rotation=0)                
                    cbar = ax.collections[0].colorbar
                    cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                    cbar.set_ticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], update_ticks=True)
                    cbar.update_ticks() 
                    plt.savefig('_post/{0}/pos_pic_file/_{1}_ms.png'.format(str(self.tgt_dirs[idx]), file_time),bbox_inches='tight')
                    #plt.show()
                    plt.close()

                else:
                    rel_pos_square = (df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)'])
                    rel_pos_unit_x = df_pos['pos_x (m)'] / np.sqrt(rel_pos_square)
                    ptcl_count = np.zeros((division_number_to_radius, division_number_to_radian)).astype(int)

                    # particle spatial distribution at each time step
                    for idx_ptcl in range (len(df_pos)):

                        r_num = np.where(rel_pos_square[idx_ptcl] <= vals_r_square)[0][0]
                        rad_num = np.where(np.arccos(np.clip(((rel_pos_unit_x[idx_ptcl] * 1.0) + (0 * 0) + (0 * 0)), -1, 1)) <= vals_rad_sum)[0][0]
                        ptcl_count[r_num][rad_num] += 1.0

                    fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(3.0,3.0), dpi=300)
                    normed_ptcl_count = ptcl_count / float(np.max(ptcl_count))
                    colors = temp_cmap(normed_ptcl_count)

                    for idx_r in range(division_number_to_radius):

                        ax.bar(x=vals_rad_base,
                            width=vals_rad_norm, bottom=vals_r_base[idx_r], height=vals_r_norm[idx_r],
                            color=colors[idx_r], edgecolor='w', linewidth=0.0, align="edge")

                    ax.set_xlim([0, np.pi])
                    ax.set_xticks([0, np.pi*1/4, np.pi*2/4, np.pi*3/4, np.pi*4/4])
                    ax.set_xticklabels([0, r'$\frac{1}{4}$'+'$\it{\u03c0}$', r'$\frac{1}{2}$'+'$\it{\u03c0}$', r'$\frac{3}{4}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
                    ax.text(-0.1, 0.4, 'Radian  $\it{\u03b8}$', rotation=90, transform=ax.transAxes)
                    ax.set_ylim([0, 1.0])
                    ax.set_yticks([0, 0.5, 1.0])
                    ax.set_yticklabels([0, 12.5, 25])
                    ax.text(0.5, 0.05, 'Radius  $\it{mm}$', rotation=0, transform=ax.transAxes)
                    ax.grid(False)
                    ax.tick_params(which='both', left='on', top='on', direction='in')

                    normed_data = Normalize(vmin=np.min(normed_ptcl_count), vmax=np.max(normed_ptcl_count))
                    mappable = ScalarMappable(cmap=temp_cmap, norm=normed_data)
                    mappable._A = []
                    cbar = plt.colorbar(mappable, ax=ax, pad=0.12, shrink=0.5)
                    cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                    cbar.ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                    cbar.outline.set_visible(False)

                    plt.box(on=None)
                    plt.savefig('_post/{0}/pos_pic_file/_{1}_ms.png'.format(str(self.tgt_dirs[idx]), file_time),bbox_inches='tight')
                    #plt.show()
                    plt.close()


    def cell_layer_pos_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            pos_measurement = int(self.measurement[idx] * 1000)                
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), pos_measurement))   
            df_pos = pd.DataFrame()  
            half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
            ptcl_num_in_each_cell = np.zeros(half_cubic_length)
            sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
            volume_in_each_cell = np.diff(4.0 / 3.0 * np.pi * np.power(np.array(range(0,sphere_radius+1)),3))
            cell_number = np.append(np.array([2**3]), np.diff(np.power(np.arange(2, 50+2, 2), 3)))
            ptcl_num_in_each_cell = np.zeros(sphere_radius)
        
            
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
                
                # average velocity in each cell layer at each tine step
                for idx_ptcl in range (len(df_pos)):

                    temp_pos = [df_pos.at[idx_ptcl, 'pos_x (mm)'], df_pos.at[idx_ptcl, 'pos_y (mm)'], df_pos.at[idx_ptcl, 'pos_z (mm)']]
                    dis_x = (half_cubic_length - 1 - temp_pos[0]).astype(int) if (temp_pos[0]<half_cubic_length) else (temp_pos[0] - half_cubic_length).astype(int)
                    dis_y = (half_cubic_length - 1 - temp_pos[1]).astype(int) if (temp_pos[1]<half_cubic_length) else (temp_pos[1] - half_cubic_length).astype(int)
                    dis_z = (half_cubic_length - 1 - temp_pos[2]).astype(int) if (temp_pos[2]<half_cubic_length) else (temp_pos[2] - half_cubic_length).astype(int)
                    in_cell = max([dis_x, dis_y, dis_z])
                    ptcl_num_in_each_cell[in_cell] += 1        

            else:        
                rel_pos = ((np.sqrt((df_pos['pos_x (m)'] * df_pos['pos_x (m)']) + (df_pos['pos_y (m)'] * df_pos['pos_y (m)']) + (df_pos['pos_z (m)'] * df_pos['pos_z (m)']))) * 1000).astype(int)
                for idx_ptcl in range (len(df_pos)):
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
            else:
                output_ave_pos_in_each_cell = np.divide(ptcl_num_in_each_cell, volume_in_each_cell, out=np.zeros_like(ptcl_num_in_each_cell), where=ptcl_num_in_each_cell!=0)


            plt.figure(figsize=(4.5,3.0),dpi=300)
            if(self.cubic_cell_bool[idx]):
                plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('cell layer')
            else:
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell/np.sum(output_ave_pos_in_each_cell), c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('layer')
            plt.ylabel('probability')
            plt.xlim(xmin=5,xmax=25)
            plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
            plt.ylim(ymin=0,ymax=0.1)
            plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], [0, 0.02, 0.04, 0.06, 0.08, 0.1],rotation=0)
            plt.savefig('_post/{0}/{1}_pos_cell_{2}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), pos_measurement), bbox_inches='tight')
            #plt.show()
            plt.close()


    def accumulated_cell_layer_pos_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            count_1 = 0
            count_2 = 0
            half_cubic_length = int((self.cell_length[idx] * 1000 * 0.5) + 0.5)
            cell_number = np.append(np.array([2**3]), np.diff(np.power(np.arange(2, 50+2, 2), 3)))
            sphere_radius = int((self.cell_length[idx] * 0.5 * 1000) + 0.5)
            volume_in_each_cell = np.diff(4.0 / 3.0 * np.pi * np.power(np.array(range(0,sphere_radius+1)),3))

            if(self.cubic_cell_bool[idx]):
                ptcl_num_in_each_cell_1 = np.zeros(half_cubic_length)
                ptcl_num_in_each_cell_2 = np.zeros(half_cubic_length)
            else:
                ptcl_num_in_each_cell_1 = np.zeros(sphere_radius)
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
                                    ptcl_num_in_each_cell_1[in_cell] += 1  
                                if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                    ptcl_num_in_each_cell_2[in_cell] += 1  
                                    
                            else:

                                if ((int(file_time) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                    ptcl_num_in_each_cell_1[rel_pos[idx_ptcl]] += 1
                                if ((int(file_time) + int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) % (int(self.turn_on_for_accum[idx] * 1000) + int(self.turn_off_for_accum[idx] * 1000)) == 0):
                                    ptcl_num_in_each_cell_2[rel_pos[idx_ptcl]] += 1

                        
            for idx_cell in range (4):
                ptcl_num_in_each_cell_1[4] += ptcl_num_in_each_cell_1[idx_cell]
                ptcl_num_in_each_cell_2[4] += ptcl_num_in_each_cell_2[idx_cell]
                cell_number[4] += cell_number[idx_cell]
                volume_in_each_cell[4] += volume_in_each_cell[idx_cell]
                ptcl_num_in_each_cell_1[idx_cell] = 0
                ptcl_num_in_each_cell_2[idx_cell] = 0
                cell_number[idx_cell] = 0
                volume_in_each_cell[idx_cell] = 0


            if(self.cubic_cell_bool[idx]):
                output_ave_pos_in_each_cell_1 = np.divide(ptcl_num_in_each_cell_1, cell_number, out=np.zeros_like(ptcl_num_in_each_cell_1), where=ptcl_num_in_each_cell_1!=0)
                output_ave_pos_in_each_cell_2 = np.divide(ptcl_num_in_each_cell_2, cell_number, out=np.zeros_like(ptcl_num_in_each_cell_2), where=ptcl_num_in_each_cell_2!=0)
            else: 
                output_ave_pos_in_each_cell_1 = np.divide(ptcl_num_in_each_cell_1, volume_in_each_cell, out=np.zeros_like(ptcl_num_in_each_cell_1), where=ptcl_num_in_each_cell_1!=0)
                output_ave_pos_in_each_cell_2 = np.divide(ptcl_num_in_each_cell_2, volume_in_each_cell, out=np.zeros_like(ptcl_num_in_each_cell_2), where=ptcl_num_in_each_cell_2!=0)


            plt.figure(figsize=(4.5,3),dpi=300)
            if(self.cubic_cell_bool[idx]):
                plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell_1/np.sum(output_ave_pos_in_each_cell_1), c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('cell layer')
            else:
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell_1/np.sum(output_ave_pos_in_each_cell_1), c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('layer')   

            plt.ylabel('probability')
            plt.xlim(xmin=5,xmax=25)
            plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
            plt.ylim(ymin=0,ymax=0.1)
            plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], [0, 0.02, 0.04, 0.06, 0.08, 0.1],rotation=0)        
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_excitation_pos_cell.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_1, count_1 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()


            plt.figure(figsize=(4.5,3),dpi=300)
            if(self.cubic_cell_bool[idx]):
                plt.plot(np.array(range(1,half_cubic_length + 1)), output_ave_pos_in_each_cell_2/np.sum(output_ave_pos_in_each_cell_2), c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('cell layer')
            else:
                plt.plot(np.array(range(1,sphere_radius+1)), output_ave_pos_in_each_cell_2/np.sum(output_ave_pos_in_each_cell_2), c=temp_cmap(1.0), ls='-', lw = 1.25)
                plt.xlabel('layer')    

            plt.ylabel('probability')
            plt.xlim(xmin=5,xmax=25)
            plt.xticks([5, 10, 15, 20, 25], [5, 10, 15, 20, 25],rotation=0)
            plt.ylim(ymin=0,ymax=0.1)
            plt.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], [0, 0.02, 0.04, 0.06, 0.08, 0.1],rotation=0)        
            plt.savefig('_post/{0}/{1}_{2}_accum_{3}_ptcls_end_cooling_pos_cell.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), count_2, count_2 * len(df)), bbox_inches='tight')
            #plt.show()
            plt.close()

            
    def count_collisions(self):

        for idx in range(len(self.tgt_dirs)):

            csv_files_path = path.Path(cr_path / self.tgt_dirs[idx] / 'OUTPUT_2')
            csv_files = np.array(list(csv_files_path.glob('_*_ms.csv')))
            coll_ptcl = np.zeros(len(csv_files)) 
            coll_wall = np.zeros(len(csv_files))
            elapsed_time = np.zeros(len(csv_files)) 
            idx_csv_start_cooling = 0
            idx_csv_start_cooling_Flag = True 

            # search and process each csv file
            for idx_csv, val_csv in enumerate (csv_files):

                # read csv file
                df = pd.read_csv(val_csv)
                coll_ptcl[idx_csv] = int(np.sum(df['Collisions_ptcl'])) if int(np.sum(df['Collisions_ptcl'])) == 0 else int(np.sum(df['Collisions_ptcl'])) / self.ptcl_num[idx]
                coll_wall[idx_csv] = int(np.sum(df['Collisions_wall'])) if int(np.sum(df['Collisions_wall'])) == 0 else int(np.sum(df['Collisions_wall'])) / self.ptcl_num[idx]
                file_time = re.findall(r'\d+', str(val_csv.name))[0]
                elapsed_time[idx_csv] = int(file_time)
                
                if (int(file_time) >= int(self.end_excitation[idx]*1000)) and (idx_csv_start_cooling_Flag == True):
                    idx_csv_start_cooling = idx_csv
                    idx_csv_start_cooling_Flag = False
                    
            coll_ptcl_cooling = coll_ptcl[idx_csv_start_cooling:]
            coll_wall_cooling = coll_wall[idx_csv_start_cooling:]
            elapsed_time_cooling = elapsed_time[idx_csv_start_cooling:]

            # showing particle average-velocity
            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(elapsed_time / 1000, coll_ptcl, c=temp_cmap(1.0), ls='solid', lw = 0.7, label='with particle')
            plt.plot(elapsed_time / 1000, coll_wall, c=temp_cmap(0.5), ls='solid', lw = 0.7, label='with boundary')
            plt.xlabel('time  $\it{s}$')
            plt.ylabel('accumulated collision times  (per particle)')        
            plt.legend(bbox_to_anchor=(1, 1.01), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_collision_accumulated.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx])),bbox_inches='tight')
            #plt.show()
            plt.close()

            # showing particle average-velocity
            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(elapsed_time / 1000, np.append(np.array([0]), np.diff(coll_ptcl)), c=temp_cmap(1.0), ls='solid', lw = 0.7, label='with particle')
            plt.plot(elapsed_time / 1000, np.append(np.array([0]), np.diff(coll_wall)), c=temp_cmap(0.5), ls='solid', lw = 0.7, label='with boundary')
            plt.xlabel('time  $\it{s}$')
            plt.ylabel('collision frequency  (per particle)')        
            plt.legend(bbox_to_anchor=(1, 1.01), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_collision_frequency_1.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx])),bbox_inches='tight')
            #plt.show()
            plt.close()

            # showing particle average-velocity
            plt.figure(figsize=(4.5,3),dpi=300)
            plt.plot(elapsed_time_cooling / 1000, np.append(np.array([0]), np.diff(coll_ptcl_cooling)), c=temp_cmap(1.0), ls='solid', lw = 0.7, label='with particle')
            plt.plot(elapsed_time_cooling / 1000, np.append(np.array([0]), np.diff(coll_wall_cooling)), c=temp_cmap(0.5), ls='solid', lw = 0.7, label='with boundary')
            plt.xlabel('time  $\it{s}$')
            plt.ylabel('collision frequency  (per particle)')  
                
            plt.xlim(xmin=5, xmax=10)
            plt.xticks([5, 6, 7, 8, 9, 10], [5, 6, 7, 8, 9, 10], rotation=0)     
            """  
            plt.ylim(ymin=0,ymax=0.5)
            plt.yticks([0,0.1,0.2,0.3,0.4,0.5], [0,0.1,0.2,0.3,0.4,0.5],rotation=0)
            plt.grid(which='minor')
            """
            plt.legend(bbox_to_anchor=(1, 1.01), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_collision_frequency_2.png'.format(str(self.tgt_dirs[idx]),str(self.output_name[idx])),bbox_inches='tight')
            #plt.show()
            plt.close()

    

