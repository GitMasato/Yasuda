from . import functions as func
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib as path
import seaborn as sns

temp_cmap = plt.get_cmap("Blues")
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams["mathtext.fontset"] = 'dejavusans'
plt.rcParams['font.size'] = 10
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
cr_path = path.Path()


class Mag_data_output:

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


    def mag_force_dist(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            mag_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), mag_measurement))
            
            x_for_with= np.logspace(np.log10(df['Mag_with_ab (m/s2)'].values.min()), np.log10(df['Mag_with_ab (m/s2)'].values.max()), self.normal_bins_size[idx]+1)
            x_for_without = np.logspace(np.log10(df['Mag_without_ab (m/s2)'].values.min()), np.log10(df['Mag_without_ab (m/s2)'].values.max()), self.normal_bins_size[idx]+1)
            plt.close()
            plt.figure(figsize=(4.5,3),dpi=300)
            
            plt.hist(df['Mag_with_ab (m/s2)'].values, bins=x_for_with, density=True, log=True, alpha=0.5, color=temp_cmap(1.0), label='with interaction')
            plt.hist(df['Mag_without_ab (m/s2)'].values, bins=x_for_without, density=True, log=True, alpha=0.7, color=temp_cmap(0.5), label='without interaction')
            plt.ylabel('probability density')
            plt.xscale('log')
            plt.xlabel('magnetic force  $\it{m/s}$$^2$')
            #plt.xlim(xmin=10**-7,xmax=10**-2)
            #plt.ylim(ymin=10**-1,ymax=10**3)
            plt.legend(bbox_to_anchor=(1, 1.01), loc='lower right', borderaxespad=0)
            plt.savefig('_post/{0}/{1}_mag_force_dist_prob_{2:0>5}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), mag_measurement), bbox_inches='tight')
            plt.close()


    def mag_dir_diff_1(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            if(self.cubic_cell_bool[idx]):

                mag_measurement = int(self.measurement[idx] * 1000)
                df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), mag_measurement)) 
                
                pos_x = df['pos_x (m)'].values - 0.005
                pos_y = df['pos_y (m)'].values - 0.005
                pos_z = df['pos_z (m)'].values - 0.005
                mag_x_unit = df['Mag_without_x (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
                mag_y_unit = df['Mag_without_y (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
                mag_z_unit = df['Mag_without_z (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
                base_vector_x = (mag_x_unit * 0.0) + 0.
                base_vector_y = (mag_z_unit / np.sqrt((mag_z_unit * mag_z_unit) + (mag_y_unit * mag_y_unit))) + 0.
                base_vector_z = (- mag_y_unit / np.sqrt((mag_z_unit * mag_z_unit) + (mag_y_unit * mag_y_unit))) + 0.
                mag_x_unit_with_interaction = df['Mag_with_x (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
                mag_y_unit_with_interaction = df['Mag_with_y (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
                mag_z_unit_with_interaction = df['Mag_with_z (m/s2)'].values / df['Mag_with_ab (m/s2)'].values  
                
                inner_dot_with_without = (mag_x_unit_with_interaction * mag_x_unit) + (mag_y_unit_with_interaction * mag_y_unit) + (mag_z_unit_with_interaction * mag_z_unit)
                mag_angle_diff_theta = np.arccos(np.clip(inner_dot_with_without, -1, 1)) / np.pi
                projected_mag_x_with_interaction = mag_x_unit_with_interaction - (inner_dot_with_without * mag_x_unit)
                projected_mag_y_with_interaction = mag_y_unit_with_interaction - (inner_dot_with_without * mag_y_unit)
                projected_mag_z_with_interaction = mag_z_unit_with_interaction - (inner_dot_with_without * mag_z_unit)
                projected_mag_x_unit_with_interaction = projected_mag_x_with_interaction / np.sqrt((projected_mag_x_with_interaction * projected_mag_x_with_interaction) + (projected_mag_y_with_interaction * projected_mag_y_with_interaction) + (projected_mag_z_with_interaction * projected_mag_z_with_interaction))
                projected_mag_y_unit_with_interaction = projected_mag_y_with_interaction / np.sqrt((projected_mag_x_with_interaction * projected_mag_x_with_interaction) + (projected_mag_y_with_interaction * projected_mag_y_with_interaction) + (projected_mag_z_with_interaction * projected_mag_z_with_interaction))
                projected_mag_z_unit_with_interaction = projected_mag_z_with_interaction / np.sqrt((projected_mag_x_with_interaction * projected_mag_x_with_interaction) + (projected_mag_y_with_interaction * projected_mag_y_with_interaction) + (projected_mag_z_with_interaction * projected_mag_z_with_interaction))
                mag_angle_diff_phi = np.zeros(len(df))
        
                temp_ang = np.radians(np.array(range(0,360)) * 1.0)
                dia_line_pos = np.array([0.0, self.cell_length[idx], 0.0])
                dia_line_vector_unit = np.array([(1./np.sqrt(3.)), - (1./np.sqrt(3.)), (1.0/np.sqrt(3.))])
                closest_vec_x = np.zeros(len(df))
                closest_vec_y = np.zeros(len(df))
                closest_vec_z = np.zeros(len(df))
                epsilon = 0.0001
                ptcl_not_contact = []

                for idx_ptcl in range(len(df)):

                    ptcl_pos = np.array([pos_x[idx_ptcl], pos_y[idx_ptcl], pos_z[idx_ptcl]])
                    r_pos = dia_line_pos - ptcl_pos
                    closest_tau = 0.
                    closest_distance = 0.
                    closest_temp_vec = np.array([0., 0., 0.])
                    first_contact_frag = True
                    first_non_contact_frag = True

                    for idx_ang, val_ang in enumerate (temp_ang):
                        
                        ROT = func.rot_mat(val_ang, [mag_x_unit[idx_ptcl], mag_y_unit[idx_ptcl], mag_z_unit[idx_ptcl]])
                        temp_vector = np.dot(ROT, [[base_vector_x[idx_ptcl]], [base_vector_y[idx_ptcl]], [base_vector_z[idx_ptcl]]])
                        ptcl_line_vector_unit = np.array([temp_vector[0][0], temp_vector[1][0], temp_vector[2][0]])
                        inner_dot_line_vectors = np.dot(dia_line_vector_unit, ptcl_line_vector_unit)
                        sigma = np.dot(r_pos, (dia_line_vector_unit - (inner_dot_line_vectors*ptcl_line_vector_unit))) / ((inner_dot_line_vectors*inner_dot_line_vectors) - 1.)
                        tau = np.dot(r_pos, ((inner_dot_line_vectors*dia_line_vector_unit) - ptcl_line_vector_unit)) / ((inner_dot_line_vectors*inner_dot_line_vectors) - 1.)                
                        distance_lines = dia_line_pos + (sigma * dia_line_vector_unit) - ptcl_pos - (tau * ptcl_line_vector_unit) 
                        distance = np.linalg.norm(distance_lines)
                        
                        if (tau >= 0):
                            if(distance < epsilon):                        
                                if (first_contact_frag == True):

                                    closest_vec_x[idx_ptcl] = ptcl_line_vector_unit[0]
                                    closest_vec_y[idx_ptcl] = ptcl_line_vector_unit[1]
                                    closest_vec_z[idx_ptcl] = ptcl_line_vector_unit[2]
                                    closest_tau = tau
                                    first_contact_frag = False
                                
                                else:
                                    if (tau < closest_tau):
                                        
                                        closest_vec_x[idx_ptcl] = ptcl_line_vector_unit[0]
                                        closest_vec_y[idx_ptcl] = ptcl_line_vector_unit[1]
                                        closest_vec_z[idx_ptcl] = ptcl_line_vector_unit[2]
                                        closest_tau = tau

                            else:
                                if (first_non_contact_frag == True):

                                    closest_temp_vec[0] = ptcl_line_vector_unit[0]
                                    closest_temp_vec[1] = ptcl_line_vector_unit[1]
                                    closest_temp_vec[2] = ptcl_line_vector_unit[2]
                                    closest_distance = distance
                                    first_non_contact_frag = False

                                else:
                                    if (distance < closest_distance):

                                        closest_temp_vec[0] = ptcl_line_vector_unit[0]
                                        closest_temp_vec[1] = ptcl_line_vector_unit[1]
                                        closest_temp_vec[2] = ptcl_line_vector_unit[2]
                                        closest_distance = distance

                    if (first_contact_frag == True):
                        
                        closest_vec_x[idx_ptcl] = closest_temp_vec[0]
                        closest_vec_y[idx_ptcl] = closest_temp_vec[1]
                        closest_vec_z[idx_ptcl] = closest_temp_vec[2] 
                        ptcl_not_contact.append(idx_ptcl)


                for idx_ptcl in range(len(df)):

                    vector_A = np.array([closest_vec_x[idx_ptcl], closest_vec_y[idx_ptcl], closest_vec_z[idx_ptcl]])
                    vector_B = np.array([projected_mag_x_unit_with_interaction[idx_ptcl], projected_mag_y_unit_with_interaction[idx_ptcl], projected_mag_z_unit_with_interaction[idx_ptcl]])
                    vector_C = np.array([mag_x_unit[idx_ptcl], mag_y_unit[idx_ptcl], mag_z_unit[idx_ptcl]])
                    inner_dot = np.dot(vector_A,vector_B)
                    outer_cross = np.cross(vector_A,vector_B) 
                    mag_angle_diff_phi[idx_ptcl] = np.where(((np.dot(outer_cross, vector_C) / np.linalg.norm(outer_cross)) > 0.9), np.arccos(np.clip(inner_dot, -1, 1)) / np.pi, - np.arccos(np.clip(inner_dot, -1, 1)) / np.pi)
        

                temp_bins = (np.arccos(1. - np.linspace(0., 2., self.normal_bins_size[idx]+1))) / np.pi
                weights_theta = plt.hist(mag_angle_diff_theta, bins=temp_bins, density=False)
                weights_phi = plt.hist(mag_angle_diff_phi, bins=self.normal_bins_size[idx], density=False)        
                plt.close()

                h = sns.JointGrid(mag_angle_diff_theta, mag_angle_diff_phi, height=5, ratio=2, space=0.3)
                h.ax_joint.cla() 
                h.ax_marg_x.cla()
                h.ax_marg_y.cla()
                h.ax_joint.scatter(mag_angle_diff_theta, mag_angle_diff_phi, s=1.0, cmap=temp_cmap(0.5))

                h.ax_joint.set_xlim([0, 1.0])
                h.ax_joint.set_xticks([0, 0.5, 1.0])
                h.ax_joint.set_xticklabels([0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
                h.ax_joint.set_ylim([-1.0, 1.0])
                h.ax_joint.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
                h.ax_joint.set_yticklabels(['-$\it{\u03c0}$', '-'+r'$\frac{1}{2}$'+'$\it{\u03c0}$', 0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
                h.set_axis_labels('angle $\it{\u03b8}$$\it{_{ res}}$', 'angle $\it{\u03c6}$$\it{_{ res}}$', fontsize=10)
                h.fig.set_figwidth(4.0)
                h.fig.set_figheight(4.0)
                
                ax1 = h.ax_marg_x
                ax2 = h.ax_marg_y
                ax1.hist(mag_angle_diff_theta, bins=temp_bins, density=False, weights=np.ones(len(mag_angle_diff_theta))/np.sum(weights_theta[0]), color=temp_cmap(0.5))
                ax2.hist(mag_angle_diff_phi, bins=self.normal_bins_size[idx], orientation="horizontal", density=False, weights=np.ones(len(mag_angle_diff_phi))/np.sum(weights_phi[0]), color=temp_cmap(0.5))        

                #ax1.set_ylim(0, 0.02)
                #ax1.yaxis.set_ticks([0, 0.01, 0.02])
                ax1.set_ylabel('probability\ndensity')
                plt.setp(ax1.get_xticklabels(), visible=False)        
                #ax2.set_xlim(0, 0.02)
                #ax2.xaxis.set_ticks([0, 0.01, 0.02])
                ax2.set_xlabel('probability\ndensity')
                plt.setp(ax2.get_yticklabels(), visible=False)        
                #plt.title('angle\ndifference', fontdict=dict(fontsize=10), position=(0.35, 1.3), bbox=dict(facecolor='w', alpha=0.1))
                h.savefig('_post/{0}/{1}_mag_force_angle_diff_prob_1_{2:0>5}_ms_non_contact_{3}.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), mag_measurement, len(ptcl_not_contact)), bbox_inches='tight', dpi=300)
                plt.close()


    def mag_dir_diff_2(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            mag_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), mag_measurement)) 
            
            mag_x_unit = df['Mag_without_x (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_y_unit = df['Mag_without_y (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_z_unit = df['Mag_without_z (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_x_unit_with_interaction = df['Mag_with_x (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_y_unit_with_interaction = df['Mag_with_y (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_z_unit_with_interaction = df['Mag_with_z (m/s2)'].values / df['Mag_with_ab (m/s2)'].values  

            inner_dot_with_without = (mag_x_unit_with_interaction * mag_x_unit) + (mag_y_unit_with_interaction * mag_y_unit) + (mag_z_unit_with_interaction * mag_z_unit)
            mag_angle_diff_theta = np.arccos(np.clip(inner_dot_with_without, -1, 1)) / np.pi
            temp_bins = (np.arccos(1. - np.linspace(0., 2., self.normal_bins_size[idx]+1))) / np.pi
            weights_theta = plt.hist(mag_angle_diff_theta, bins=temp_bins, density=False)
            plt.close()

            plt.figure(figsize=(4.0,4.0),dpi=300)
            plt.hist(mag_angle_diff_theta, bins=temp_bins, density=False, weights=np.ones(len(mag_angle_diff_theta))/np.sum(weights_theta[0]), color=temp_cmap(0.5))
            plt.xlabel('angle $\it{\u03b8}$$_{ res}$')
            plt.ylabel('probability density')
            plt.xlim(xmin=0,xmax=1.0)
            plt.xticks([0, 0.5, 1.0], [0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'],rotation=0)
            #plt.ylim(ymin=0,ymax=0.04)
            #plt.yticks([0, 0.01, 0.02, 0.03, 0.04], [0, 0.01, 0.02, 0.03, 0.04],rotation=0)
            
            plt.savefig('_post/{0}/{1}_mag_force_angle_diff_prob_3_{2:0>5}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), mag_measurement), bbox_inches='tight')
            #plt.show()
            plt.close()


    def mag_dir_dist_1(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):

            mag_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), mag_measurement))
            
            mag_x_unit_with = df['Mag_with_x (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_y_unit_with = df['Mag_with_y (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_z_unit_with = df['Mag_with_z (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_phi_with = np.arctan2(mag_y_unit_with, mag_x_unit_with) / np.pi
            mag_theta_with = np.arccos(np.clip(mag_z_unit_with, -1, 1)) / np.pi

            temp_bins = (np.arccos(1. - np.linspace(0., 2., self.normal_bins_size[idx]+1))) / np.pi
            weights_theta_with = plt.hist(mag_theta_with, bins=temp_bins, density=False)
            weights_phi_with = plt.hist(mag_phi_with, bins=self.normal_bins_size[idx], density=False)        
            plt.close()


            h = sns.JointGrid(mag_theta_with, mag_phi_with, height=5, ratio=2, space=0.3)
            h.ax_joint.cla() 
            h.ax_marg_x.cla()
            h.ax_marg_y.cla()
            h.ax_joint.scatter(mag_theta_with, mag_phi_with, s=1.0, cmap=temp_cmap(0.5))

            h.ax_joint.set_xlim([0, 1.0])
            h.ax_joint.set_xticks([0, 0.5, 1.0])
            h.ax_joint.set_xticklabels([0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
            h.ax_joint.set_ylim([-1.0, 1.0])
            h.ax_joint.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
            h.ax_joint.set_yticklabels(['-$\it{\u03c0}$', '-'+r'$\frac{1}{2}$'+'$\it{\u03c0}$', 0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
            h.set_axis_labels('angle $\it{\u03b8}$$\it{_{ z}}$', 'angle $\it{\u03c6}$$\it{_{ x-y}}$', fontsize=10)
            h.fig.set_figwidth(4.0)
            h.fig.set_figheight(4.0)
            
            ax1 = h.ax_marg_x
            ax2 = h.ax_marg_y
            ax1.hist(mag_theta_with, bins=temp_bins, density=False, weights=np.ones(len(mag_theta_with))/np.sum(weights_theta_with[0]), color=temp_cmap(0.5))
            ax2.hist(mag_phi_with, bins=self.normal_bins_size[idx], orientation="horizontal", density=False, weights=np.ones(len(mag_phi_with))/np.sum(weights_phi_with[0]), color=temp_cmap(0.5))        

            #ax1.set_ylim(0, 0.02)
            #ax1.yaxis.set_ticks([0, 0.01, 0.02])
            ax1.set_ylabel('probability\ndensity')
            plt.setp(ax1.get_xticklabels(), visible=False)        
            #ax2.set_xlim(0, 0.02)
            #ax2.xaxis.set_ticks([0, 0.01, 0.02])
            ax2.set_xlabel('probability\ndensity')
            plt.setp(ax2.get_yticklabels(), visible=False)        
            #plt.title('with\ninteraction', fontdict=dict(fontsize=10), position=(0.35, 1.3), bbox=dict(facecolor='w', alpha=0.1))
            h.savefig('_post/{0}/{1}_mag_force_dir_dist_prob_1_with_interaction_{2:0>5}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), mag_measurement), bbox_inches='tight', dpi=300)
            plt.close()


            mag_x_unit_without = df['Mag_without_x (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_y_unit_without = df['Mag_without_y (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_z_unit_without = df['Mag_without_z (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_phi_without = np.arctan2(mag_y_unit_without, mag_x_unit_without) / np.pi
            mag_theta_without = np.arccos(np.clip(mag_z_unit_without, -1, 1)) / np.pi

            temp_bins = (np.arccos(1. - np.linspace(0., 2., self.normal_bins_size[idx]+1))) / np.pi
            weights_theta_without = plt.hist(mag_theta_without, bins=temp_bins, density=False)
            weights_phi_without = plt.hist(mag_phi_without, bins=self.normal_bins_size[idx], density=False)        
            plt.close()
            
            h = sns.JointGrid(mag_theta_without, mag_phi_without, height=5, ratio=2, space=0.3)
            h.ax_joint.cla() 
            h.ax_marg_x.cla()
            h.ax_marg_y.cla()
            h.ax_joint.scatter(mag_theta_without, mag_phi_without, s=1.0, cmap=temp_cmap(0.5))

            h.ax_joint.set_xlim([0, 1.0])
            h.ax_joint.set_xticks([0, 0.5, 1.0])
            h.ax_joint.set_xticklabels([0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
            h.ax_joint.set_ylim([-1.0, 1.0])
            h.ax_joint.set_yticks([-1.0, -0.5, 0, 0.5, 1.0])
            h.ax_joint.set_yticklabels(['-$\it{\u03c0}$', '-'+r'$\frac{1}{2}$'+'$\it{\u03c0}$', 0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'])
            h.set_axis_labels('angle $\it{\u03b8}$$\it{_{ z}}$', 'angle $\it{\u03c6}$$\it{_{ x-y}}$', fontsize=10)
            h.fig.set_figwidth(4.0)
            h.fig.set_figheight(4.0)
            
            ax1 = h.ax_marg_x
            ax2 = h.ax_marg_y
            ax1.hist(mag_theta_without, bins=temp_bins, density=False, weights=np.ones(len(mag_theta_without))/np.sum(weights_theta_without[0]), color=temp_cmap(0.5))
            ax2.hist(mag_phi_without, bins=self.normal_bins_size[idx], orientation="horizontal", density=False, weights=np.ones(len(mag_phi_without))/np.sum(weights_phi_without[0]), color=temp_cmap(0.5))        

            #ax1.set_ylim(0, 0.02)
            #ax1.yaxis.set_ticks([0, 0.01, 0.02])
            ax1.set_ylabel('probability\ndensity')
            plt.setp(ax1.get_xticklabels(), visible=False)        
            #ax2.set_xlim(0, 0.02)
            #ax2.xaxis.set_ticks([0, 0.01, 0.02])
            ax2.set_xlabel('probability\ndensity')
            plt.setp(ax2.get_yticklabels(), visible=False)        
            #plt.title('without\ninteraction', fontdict=dict(fontsize=10), position=(0.35, 1.3), bbox=dict(facecolor='w', alpha=0.1))
            h.savefig('_post/{0}/{1}_mag_force_dir_dist_prob_1_without_interaction_{2:0>5}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), mag_measurement), bbox_inches='tight', dpi=300)
            plt.close()


    def mag_dir_dist_2(self):

        # search and process in each directory
        for idx in range(len(self.tgt_dirs)):
        
            mag_measurement = int(self.measurement[idx] * 1000)
            df = pd.read_csv('{0}/OUTPUT_2/_{1:0>5}_ms.csv'.format(str(self.tgt_dirs[idx]), mag_measurement))
            
            mag_x_unit_with = df['Mag_with_x (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_y_unit_with = df['Mag_with_y (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_z_unit_with = df['Mag_with_z (m/s2)'].values / df['Mag_with_ab (m/s2)'].values
            mag_phi_with = np.arctan2(mag_y_unit_with, mag_x_unit_with) / np.pi
            mag_theta_with = np.arccos(np.clip(mag_z_unit_with, -1, 1)) / np.pi

            plt.figure(figsize=(4.0,4.0),dpi=300)
            plt.scatter(mag_theta_with, mag_phi_with, s=1, cmap=temp_cmap(0.5))
            plt.xlabel('angle $\it{\u03b8}$$\it{_{ z}}$')
            plt.ylabel('angle $\it{\u03c6}$$\it{_{ x-y}}$')
            plt.xlim(xmin=0,xmax=1.0)
            plt.xticks([0, 0.5, 1.0], [0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'],rotation=0)
            plt.ylim(ymin=-1.0,ymax=1.0)
            plt.yticks([-1.0, -0.5, 0, 0.5, 1.0], ['-$\it{\u03c0}$', '-'+r'$\frac{1}{2}$'+'$\it{\u03c0}$', 0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'],rotation=0)
            plt.savefig('_post/{0}/{1}_mag_force_dir_dist_prob_2_with_interaction_{2:0>5}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), mag_measurement), bbox_inches='tight')
            #plt.show()
            plt.close()


            mag_x_unit_without = df['Mag_without_x (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_y_unit_without = df['Mag_without_y (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_z_unit_without = df['Mag_without_z (m/s2)'].values / df['Mag_without_ab (m/s2)'].values
            mag_phi_without = np.arctan2(mag_y_unit_without, mag_x_unit_without) / np.pi
            mag_theta_without = np.arccos(np.clip(mag_z_unit_without, -1, 1)) / np.pi

            plt.figure(figsize=(4.0,4.0),dpi=300)
            plt.scatter(mag_theta_without, mag_phi_without, s=1, cmap=temp_cmap(0.5))
            plt.xlabel('angle $\it{\u03b8}$$\it{_{ z}}$')
            plt.ylabel('angle $\it{\u03c6}$$\it{_{ x-y}}$')
            plt.xlim(xmin=0,xmax=1.0)
            plt.xticks([0, 0.5, 1.0], [0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'],rotation=0)
            plt.ylim(ymin=-1.0,ymax=1.0)
            plt.yticks([-1.0, -0.5, 0, 0.5, 1.0], ['-$\it{\u03c0}$', '-'+r'$\frac{1}{2}$'+'$\it{\u03c0}$', 0, r'$\frac{1}{2}$'+'$\it{\u03c0}$', '$\it{\u03c0}$'],rotation=0)
            plt.savefig('_post/{0}/{1}_mag_force_dir_dist_prob_2_without_interaction_{2:0>5}_ms.png'.format(str(self.tgt_dirs[idx]), str(self.output_name[idx]), mag_measurement), bbox_inches='tight')
            #plt.show()
            plt.close()
   
        