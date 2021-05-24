from _MEGraMa_modules import MEGraMa_pos as Pos
from _MEGraMa_modules import MEGraMa_velo as Velo
from _MEGraMa_modules import MEGraMa_mag as Mag
from _MEGraMa_modules import MEGraMa_plural as Plural
import numpy as np
import pathlib as path
ser_dirs = np.array(['sph_5s_On20ms_Off80ms_02', 'sph_5s_On20ms_Off80ms_055', 'sph_5s_On20ms_Off80ms_09'], dtype=str)      # list of search keyword
'''
'box_5s_On20ms_Off80ms_02', 'box_5s_On20ms_Off80ms_055', 'box_5s_On20ms_Off80ms_09', \
'sph_5s_On20ms_Off80ms_02', 'sph_5s_On20ms_Off80ms_055', 'sph_5s_On20ms_Off80ms_09', \
'5s_On01ms_Off99ms_1A', '5s_On05ms_Off95ms_1A', '5s_On10ms_Off90ms_1A', '5s_On20ms_Off80ms_1A', \
'5s_On40ms_Off60ms_1A', '5s_On60ms_Off40ms_1A', '5s_On80ms_Off20ms_1A', '5s_On100ms_Off0ms_1A', \
'5s_On01ms_Off99ms_2A', '5s_On05ms_Off95ms_2A', '5s_On10ms_Off90ms_2A', '5s_On20ms_Off80ms_2A', \
'5s_On40ms_Off60ms_2A', '5s_On60ms_Off40ms_2A', '5s_On80ms_Off20ms_2A', '5s_On100ms_Off0ms_2A', \
'5s_On01ms_Off99ms_3A', '5s_On05ms_Off95ms_3A', '5s_On10ms_Off90ms_3A', '5s_On20ms_Off80ms_3A', \
'5s_On40ms_Off60ms_3A', '5s_On60ms_Off40ms_3A', '5s_On80ms_Off20ms_3A', '5s_On100ms_Off0ms_3A', \
'1570_ini_5s_On20ms_Off80ms', '3140_5s_On20ms_Off80ms_1', '4710_5s_On20ms_Off80ms_1', \
'10s_On05ms_Off95ms_1A', '20s_On01ms_Off99ms_2A', '45s_On01ms_Off99ms_1A' \
'''

output_ave_velo = False
output_velo_dist = False
output_accumulated_velo_dist = False
output_velo_dist_3_dir = False
output_maxwellian_velo_dist = False
output_accumulated_maxwellian_velo_dist = False
output_velo_dir = False
output_accumulated_velo_dir = False
output_cell_layer_velo_dist = False
output_accumulated_cell_layer_velo_dist = False

output_pos_dist = False
output_accumulated_pos_dist = False
output_pos_dist_movie = False
output_cell_layer_pos_dist = False
output_accumulated_cell_layer_pos_dist = False
output_count_collisions = False

output_mag_force_dist = False
output_mag_dir_diff_1 = False
output_mag_dir_diff_2 = False
output_mag_dir_dist_1 = False
output_mag_dir_dist_2 = False

output_plural_ave_velo = True
output_plural_saturation_velo = False
output_plural_cell_layer_velo_dist = False
output_plural_accumulated_cell_layer_velo_dist = False
output_plural_cell_layer_pos_dist = False
output_plural_accumulated_cell_layer_pos_dist = False

Output_name = np.array(['sph'])
Ptcl_num = np.array([1570])
Ptcl_dia = np.array([0.0016])           # m
Ptcl_coe_res = np.array([0.2, 0.55, 0.9])
End_excitation = np.array([5.0])        # s
Cell_length = np.array([0.05])          # m  sphere_cell: cell diameter
Turn_on_for_accum = np.array([0.02])
Turn_off_for_accum = np.array([0.08])
Measurement = np.array([0.01])          # s
Normal_bins_size = np.array([100])
Time_start_accum = np.array([4.0])
Cubic_cell_bool = np.array([False])     # s
Vibration_on = np.array([False])
Angle_xy = np.array([45.])              # angle
Angle_z = np.array([35.26])             # angle
Vib_amplitude = np.array([0.0025])      # m
Vib_frequency = np.array([10.])         # Hz
Cmap_array = np.array([1.0, 0.8, 0.6])       # 0.4, 0.2
Ls_array = np.array([':', '-', ':'])       # '-', '--', ':', '-.'
Lw_array = np.array([2.0, 1.6, 1.2])
Label_array = np.array(['0.2','0.55','0.9'])


if __name__ == '__main__':

    post_path = path.Path('_post')
    post_path.mkdir(exist_ok=True)
    cr_path = path.Path()
    cr_dirs = np.array([d_name for d_name in cr_path.iterdir() if (cr_path / d_name / 'OUTPUT_2').is_dir()])
    temp_dirs = []

    for ser_dir in ser_dirs:
        for tgt_dir in cr_dirs:
            if  (tgt_dir.name == ser_dir):
                temp_dirs.append(tgt_dir.name)
                post_dir_path = path.Path(post_path / str(tgt_dir.name))
                post_dir_path.mkdir(exist_ok=True)

    Tgt_dirs = np.array(temp_dirs)
    Pos_inst = Pos.Pos_data_output()
    Pos_inst.set_values(Tgt_dirs, Output_name, Ptcl_num, Ptcl_dia, Ptcl_coe_res, End_excitation, Cell_length, Turn_on_for_accum, Turn_off_for_accum, Measurement, Normal_bins_size, Time_start_accum, Cubic_cell_bool, Vibration_on, Angle_xy, Angle_z, Vib_amplitude, Vib_frequency, Cmap_array, Ls_array, Lw_array, Label_array)
    Velo_inst = Velo.Velo_data_output()
    Velo_inst.set_values(Tgt_dirs, Output_name, Ptcl_num, Ptcl_dia, Ptcl_coe_res, End_excitation, Cell_length, Turn_on_for_accum, Turn_off_for_accum, Measurement, Normal_bins_size, Time_start_accum, Cubic_cell_bool, Vibration_on, Angle_xy, Angle_z, Vib_amplitude, Vib_frequency, Cmap_array, Ls_array, Lw_array, Label_array)
    Plural_inst = Plural.Plural_data_output()
    Plural_inst.set_values(Tgt_dirs, Output_name, Ptcl_num, Ptcl_dia, Ptcl_coe_res, End_excitation, Cell_length, Turn_on_for_accum, Turn_off_for_accum, Measurement, Normal_bins_size, Time_start_accum, Cubic_cell_bool, Vibration_on, Angle_xy, Angle_z, Vib_amplitude, Vib_frequency, Cmap_array, Ls_array, Lw_array, Label_array)
    Mag_inst = Mag.Mag_data_output()
    Mag_inst.set_values(Tgt_dirs, Output_name, Ptcl_num, Ptcl_dia, Ptcl_coe_res, End_excitation, Cell_length, Turn_on_for_accum, Turn_off_for_accum, Measurement, Normal_bins_size, Time_start_accum, Cubic_cell_bool, Vibration_on, Angle_xy, Angle_z, Vib_amplitude, Vib_frequency, Cmap_array, Ls_array, Lw_array, Label_array)

    if ((output_plural_ave_velo  == True) or (output_plural_cell_layer_velo_dist == True) or (output_plural_cell_layer_pos_dist  == True) or (output_plural_accumulated_cell_layer_velo_dist) or (output_plural_accumulated_cell_layer_pos_dist)):

        plural_dir_path = path.Path(post_path / '_plural_data')
        plural_dir_path.mkdir(exist_ok=True)


    if output_ave_velo == True:
        Velo_inst.ave_velo()

    if output_velo_dist == True:
        Velo_inst.velo_dist()

    if output_accumulated_velo_dist == True:
        Velo_inst.accumulated_velo_dist()

    if output_velo_dist_3_dir == True:
        Velo_inst.velo_dist_3_dir()

    if output_maxwellian_velo_dist == True:
        Velo_inst.maxwellian_velo_dist()

    if output_accumulated_maxwellian_velo_dist == True:
        Velo_inst.accumulated_maxwellian_velo_dist()

    if output_velo_dir == True:
        Velo_inst.velo_dir()

    if output_accumulated_velo_dir == True:
        Velo_inst.accumulated_velo_dir()

    if output_cell_layer_velo_dist == True:
        Velo_inst.cell_layer_velo_dist()

    if output_accumulated_cell_layer_velo_dist == True:
        Velo_inst.accumulated_cell_layer_velo_dist()

    if output_pos_dist == True:
        Pos_inst.pos_dist()

    if output_accumulated_pos_dist == True:
        Pos_inst.accumulated_pos_dist()

    if output_pos_dist_movie == True:
        Pos_inst.pos_dist_movie()

    if output_cell_layer_pos_dist == True:
        Pos_inst.cell_layer_pos_dist()

    if output_accumulated_cell_layer_pos_dist == True:
        Pos_inst.accumulated_cell_layer_pos_dist()

    if output_count_collisions == True:
        Pos_inst.count_collisions()

    if output_mag_force_dist == True:
        Mag_inst.mag_force_dist()

    if output_mag_dir_diff_1 == True:
        Mag_inst.mag_dir_diff_1()

    if output_mag_dir_diff_2 == True:
        Mag_inst.mag_dir_diff_2()

    if output_mag_dir_dist_1 == True:
        Mag_inst.mag_dir_dist_1()

    if output_mag_dir_dist_2 == True:
        Mag_inst.mag_dir_dist_2()

    if output_plural_ave_velo == True:
        Plural_inst.plural_ave_velo()

    if output_plural_saturation_velo == True:
        Plural_inst.plural_saturation_velo()

    if output_plural_cell_layer_velo_dist == True:
        Plural_inst.plural_cell_layer_velo_dist()

    if output_plural_accumulated_cell_layer_velo_dist == True:
        Plural_inst.plural_accumulated_cell_layer_velo_dist()

    if output_plural_cell_layer_pos_dist == True:
        Plural_inst.plural_cell_layer_pos_dist()

    if output_plural_accumulated_cell_layer_pos_dist == True:
        Plural_inst.plural_accumulated_cell_layer_pos_dist()



