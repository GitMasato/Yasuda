#include "DEM.h"

namespace dem
{

	void DEM::RunDEM() noexcept
	{	
		double time = 0.0;
		int output_step = 0, ele_update_step = 0, mag_update_step = 0;

		#pragma omp parallel num_threads(openmp_size) if(is_openmp)
		{			
			// output initial condition
			auto kinetic_energy = ptcl.GetKineticEnergy();
			#pragma omp single
			{
				output.CreateTimeFile(total_step, kinetic_energy, time_file);
				const int avs_step = total_step / intvl_output;
				if (is_avs_output) avs.CreateAVS(avs_step, gravity, domain, avs_file, avs_key);
				if (is_ptcl_output) output.CreatePtclFile(time, ptcl_file, ptcl_key);
				if (is_vtk_output) vtk.CreateFile(time, vtk_file, vtk_key);
				output_step += intvl_output;
			}		

			/////////////////////////////////////////////////////////////////////
			// dem loop /////////////////////////////////////////////////////////
			for (int step = 0; step <= total_step; step++)
			{
				// update other object info
				#pragma omp single
				{
					obj.SetHStepVelo(time);
					obj.SetPostPos(time);
				}

				// update particle info
				ptcl.TimeIntegHVeloPosVerlet(time_step);
				ptcl.TimeIntegHAngVeloAngVerlet(time_step);
				CheckBoundary();

				// assign particles into virtual cell to find neighbor 
				#pragma omp single
				if (!ptcl.isPtclStaySizeChanged()) cell.ClearContainedPtcl();
				cell.SetPtclInCell(is_periodic[2], is_cell_z, cell_domain_z);
								
				if (is_Coulomb_ptcl || is_diele_dipole || is_mag_dipole)
				{
					#pragma omp single
					if (!ptcl.isPtclStaySizeChanged()) lr_cell.ClearContainedPtcl();
					lr_cell.SetPtclInCell(is_periodic[2], is_cell_z, cell_domain_z);
				}
							
				// calculation of electrostatic force
				if ((step == ele_update_step) && (is_ele_field))
				{
					#pragma omp single
					ele_field.SetTgtField(time, ele_start_stop_time, ele_on_off_time);

					CalcEleForce();

					#pragma omp single
					ele_update_step += intvl_ele;
				}
				
				// calculation of mgnetic force
				if ((step == mag_update_step) && (is_mag_field))
				{
					#pragma omp single
					mag_field.SetTgtField(time, mag_start_stop_time, mag_on_off_time);

					CalcMagForce();

					#pragma omp single
					mag_update_step += intvl_mag;
				}
				
				// clear collision order information   
				#pragma omp single
				coll.ClearCollOrder();

				// detect collision     
				if (is_adh) ptcl.ClearAdhForce();
				if (is_image_force)	ptcl.ClearImageForce();
				coll.DetectCollPtcl(is_periodic, domain, is_openmp, openmp_size);
				if (!is_image_force) coll.DetectCollObj(is_openmp, openmp_size);
				else coll.DetectCollObjImage(is_openmp, openmp_size, coef_image);
					
				// calculate collision  
				#pragma omp single
				{
					coll.SortCollInfo();
					if (is_roll_fri && is_adh) coll.CalcCollAdhRoll(hard_sphere_info);
					else if (is_roll_fri) coll.CalcCollRoll(hard_sphere_info);
					else if (is_adh) coll.CalcCollAdh(hard_sphere_info);
					else coll.CalcColl(hard_sphere_info);
				}

				// sum forces (with gravitational force)
				SumForceTorque();
				
				// update particle and object (velocity)  
				ptcl.TimeIntegVeloVerlet(time_step);
				ptcl.TimeIntegAngVeloVerlet(time_step);
				#pragma omp single
				obj.SetPostVelo(time);

				// output file                                                                            
				if (step == output_step)
				{
					kinetic_energy = ptcl.GetKineticEnergy();		
					#pragma omp single
					{
						output.WriteTimeFile(step, total_step, kinetic_energy, time_file);
						if (is_ptcl_output) output.CreatePtclFile(time, ptcl_file, ptcl_key);
						if (is_vtk_output) vtk.CreateFile(time, vtk_file, vtk_key);
						if (is_avs_output) avs.WriteAVS(time, gravity, domain, avs_file, avs_key);
						output_step += intvl_output;
					}
				}

				#pragma omp single
				time += time_step;
			}
			// end dem loop ///////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////
		}

		output.WriteTimeFileEnd(time_file);
		if (ptcl.GetPtclStaySize()) output.CreatePtclStayFile(ptcl_stay_last_file);
		if (ptcl.GetPtclOutSize()) output.CreatePtclOutFile(ptcl_out_last_file);
	}


	void DEM::CalcEleForce() noexcept
	{
		if (ele_field.isExtFieldOn())
		{
			ele_field.SetPtclFieldInfo(is_periodic);
			if (is_diele_force || is_diele_torque) ele_field.SetMoment(coef_ele_moment);
			
			if (is_diele_dipole)
			{
				if (is_Coulomb_ptcl) ptcl.ClearCoulombPtclForce();
				if (!is_Coulomb_ptcl) lr_ptcl_list.SetPtclInfo();
				else lr_ptcl_list.SetPtclInfoCalcCoulombPtcl(coef_ele);
				ele_lr.CalcEleDipoleInteract(coef_ele_moment, coef_ele);
				ele_lr.SetFieldFromOtherDipole(coef_ele);
			}

			if (is_Coulomb_ext_field) ele_field.CalcCoulombExtForce();
			if (is_diele_force) ele_field.CalcDieleForce();
			if (is_diele_torque) ele_field.CalcDieleTorque();
		}
		else
		{
			if (is_Coulomb_ext_field) ptcl.ClearCoulombExtForce();
			if (is_diele_force) ptcl.ClearDieleForce();
			if (is_diele_torque) ptcl.ClearDieleTorque();
			if (is_Coulomb_ptcl) ptcl.ClearCoulombPtclForce();
			if (is_Coulomb_ptcl) lr_ptcl_list.CalcCoulombPtcl(coef_ele);
		}
	}


	void DEM::CalcMagForce() noexcept
	{
		if (mag_field.isExtFieldOn())
		{
			mag_field.SetPtclFieldInfo(is_periodic);
			if (is_mag_force || is_mag_torque) mag_field.SetMoment(coef_mag_moment);

			if (is_mag_dipole)
			{
				if (!lr_ptcl_list.isLRPtclListUpdated()) lr_ptcl_list.SetPtclInfo();
				lr_ptcl_list.LRPtclListisNotUpdated();
				mag_lr.CalcMagDipoleInteract(coef_mag_moment, coef_mag);
				mag_lr.SetFieldFromOtherDipole(coef_mag);
			}

			if (is_mag_force) mag_field.CalcMagForce();
			if (is_mag_torque) mag_field.CalcMagTorque();
		}
		else
		{
			if (is_mag_force) ptcl.ClearMagForce();
			if (is_mag_torque) ptcl.ClearMagTorque();
		}
	}

	void DEM::CheckBoundary() noexcept
	{
		// check if particle reach boundary (periodic)
		if (is_periodic[0]) ptcl.CheckBoundaryXX(domain.x);
		if (is_periodic[1]) ptcl.CheckBoundaryYY(domain.y);
		if (is_periodic[2]) ptcl.CheckBoundaryZZ(domain.z);

		// check if particle reach boundary
		#pragma omp single
		{
			ptcl.PtclStaySizeisNotChanged();
			if (!is_periodic[0]) ptcl.CheckBoundaryX(domain.x);
			if (!is_periodic[1]) ptcl.CheckBoundaryY(domain.y);
			if (!is_periodic[2]) ptcl.CheckBoundaryZ(domain.z);

			if (ptcl.isPtclStaySizeChanged())
			{
				cell.AdaptToPtclSize();
				if (is_ele_field) ele_field.ModifyPtcl();
				if (is_mag_field) mag_field.ModifyPtcl();

				if (is_Coulomb_ptcl || is_diele_dipole || is_mag_dipole)
				{
					lr_cell.AdaptToPtclSize();
					lr_ptcl_list.ModifyPtcl();
				}
			}
		}
	}


	void DEM::SumForceTorque() noexcept
	{
		ptcl.StorePreExtForce();
		ptcl.AddGravityForce(gravity);
		if (is_Coulomb_ext_field) ptcl.AddCoulombExtForce();
		if (is_Coulomb_ptcl) ptcl.AddCoulombPtclForce();
		if (is_diele_force) ptcl.AddDieleForce();
		if (is_image_force) ptcl.AddImageForce();
		if (is_mag_force) ptcl.AddMagForce();
		if (is_adh) ptcl.AddAdhForce();

		if (is_cunningham) ptcl.AddAirDragCunningham(air_drag_info);
		else if (is_reynolds) ptcl.AddAirDragReynolds(air_drag_info);
		else ptcl.AddAirDrag(air_drag_info);


		ptcl.StorePreExtTorque();
		if (is_diele_torque) ptcl.AddDieleTorque();
		if (is_mag_torque) ptcl.AddMagTorque();

	}

	
	void DEM::AssociateObject() noexcept
	{
		std::cout << std::endl;
		std::cout << "associating same instance..." << std::endl;

		std::cout << "associate cell" << std::endl;
		cell.AssociatePtcl(&ptcl);
		cell.ResizeCellPtclTable();
		cell.AssociateObj(&obj);

		if (is_Coulomb_ptcl || is_diele_dipole || is_mag_dipole)
		{
			std::cout << "associate lr_cell" << std::endl;
			lr_cell.AssociatePtcl(&ptcl);
			lr_cell.ResizeCellPtclTable();
			lr_cell.AssociateObj(&obj);
		}

		std::cout << "associate coll" << std::endl;
		coll.AssociatePtcl(&ptcl);
		coll.AssociateObj(&obj);
		coll.AssociateCell(&cell);

		if (is_ele_field)
		{
			std::cout << "associate ele_field" << std::endl;
			ele_field.AssociatePtcl(&ptcl);
			ele_field.ResizeArraySize();
			ele_field.AssociateObj(&obj);
		}

		if (is_mag_field)
		{
			std::cout << "associate mag_field" << std::endl;
			mag_field.AssociatePtcl(&ptcl);
			mag_field.ResizeArraySize();
			mag_field.AssociateObj(&obj);
		}

		if (is_Coulomb_ptcl || is_diele_dipole || is_mag_dipole)
		{
			if (is_ele_field)
			{
				std::cout << "associate ele_lr" << std::endl;
				ele_lr.AssociatePtcl(&ptcl);
				ele_lr.AssociateObj(&obj);
				ele_lr.AssociateField(&ele_field);
				ele_lr.AssociateLRPtcl(&lr_ptcl_list);
			}

			if (is_mag_field)
			{
				std::cout << "associate mag_lr" << std::endl;
				mag_lr.AssociatePtcl(&ptcl);
				mag_lr.AssociateObj(&obj);
				mag_lr.AssociateField(&mag_field);
				mag_lr.AssociateLRPtcl(&lr_ptcl_list);
			}

			std::cout << "associate lr_ptcl_list" << std::endl;
			lr_ptcl_list.AssociatePtcl(&ptcl);
			lr_ptcl_list.ResizeArraySize();
			lr_ptcl_list.AssociateLRCell(&lr_cell);
		}

		if (is_avs_output)
		{
			std::cout << "associate avs" << std::endl;
			avs.AssociatePtcl(&ptcl);
			avs.AssociateObj(&obj);
		}

		std::cout << "associate output" << std::endl;
		output.AssociatePtcl(&ptcl);
		if (is_ele_field) output.AssociateEleField(&ele_field);
		if (is_mag_field) output.AssociateMagField(&mag_field);

		if (is_vtk_output)
		{
			std::cout << "associate vtk" << std::endl;
			vtk.AssociatePtcl(&ptcl);
			vtk.AssociateObj(&obj);
			if (is_ele_field) vtk.AssociateEleField(&ele_field);
			if (is_mag_field) vtk.AssociateMagField(&mag_field);
		}

		std::cout << "associate ptcl" << std::endl;
		ptcl.AssociateObj(&obj);
		std::cout << "finish associating same instance" << std::endl;

	}


	void DEM::CreateCell() noexcept
	{
		std::cout << std::endl;
		std::cout << "creating virtual cell, for collision detection..." << std::endl;
		cell.CreateCell(coll_cell_ratio, domain, is_periodic[2], is_cell_z, cell_domain_z);
		cell.SetCellInfo();
		cell.SetPeriodicEffect(domain, is_periodic);
		std::cout << "linking cell..." << std::endl;
		cell.SetLinkedSurrCell();
		std::cout << "finish creating virtual cell, for collision detection" << std::endl;
		std::cout << std::endl;

		if (is_Coulomb_ptcl || is_diele_dipole || is_mag_dipole)
		{
			std::cout << "creating virtual cell, for long range interaction..." << std::endl;
			lr_cell.CreateCell(lr_cell_ratio, domain, is_periodic[2], is_cell_z, cell_domain_z);
			lr_cell.SetCellInfo();
			lr_cell.SetPeriodicEffect(domain, is_periodic);
			std::cout << "linking cell..." << std::endl;
			lr_cell.SetLinkedSurrCell();
			std::cout << "finish creating virtual cell, for long range interaction" << std::endl;
			std::cout << std::endl;
		}
				
		//cell.ShowValue();
		//if (is_Coulomb_ptcl || is_diele_dipole || is_mag_dipole) lr_cell.ShowValue();
		
	}
	

	void DEM::SetPathInput(const int& argc, char** argv) noexcept
	{
		std::cout << std::endl;
		std::cout << "setting path to input files..." << std::endl;
		namespace fs = std::filesystem;

		for (int i = 1; i < argc; ++i)
		{
			const fs::path input = argv[i];
			const std::string input_name = input.string();

			if (std::regex_search(input_name, std::regex(".inidem")))
			{
				const fs::path p = (input.is_absolute())
					? input : fs::current_path() / input;
				ini_input_file = p.string();
			}

			if (std::regex_search(input_name, std::regex("ptcl.csv")))
			{
				const fs::path p = (input.is_absolute())
					? input : fs::current_path() / input;
				ptcl_input_file = p.string();
			}

			if (std::regex_search(input_name, std::regex(".magdem")))
			{
				const fs::path p = (input.is_absolute())
					? input : fs::current_path() / input;
				mag_input_file.push_back(p.string());
			}

			if (std::regex_search(input_name, std::regex(".eledem")))
			{
				const fs::path p = (input.is_absolute())
					? input : fs::current_path() / input;
				ele_input_file.push_back(p.string());
			}

			else if (std::regex_search(input_name, std::regex(".objdem")))
			{
				const fs::path p = (input.is_absolute())
					? input : fs::current_path() / input;
				dobj_input_file.push_back(p.string());
			}
		}

		sort(mag_input_file.begin(), mag_input_file.end());
		sort(ele_input_file.begin(), ele_input_file.end());
		sort(dobj_input_file.begin(), dobj_input_file.end());

		scope.ShowVariable(ini_input_file, "ini_input_file");
		scope.ShowVariable(ptcl_input_file, "ptcl_input_file");
		scope.ShowArray(mag_input_file, "mag_input_file");
		scope.ShowArray(ele_input_file, "ele_input_file");
		scope.ShowArray(dobj_input_file, "dobj_input_file");
		std::cout << "finish setting path to input files" << std::endl;
		std::cout << std::endl;
	}


	void DEM::SetPathOutput() noexcept
	{
		std::cout << std::endl;
		std::cout << "setting path to output files..." << std::endl;
		namespace fs = std::filesystem;

		// create directory to store output file
		fs::path p_cr = fs::current_path();
		fs::path p_ini = ini_input_file;
		std::string ini_filename = p_ini.filename().string();
		ini_filename.erase(ini_filename.end() - 7, ini_filename.end());

		fs::path p_cr_plus_ini_name = p_cr;
		p_cr_plus_ini_name.append(ini_filename);
		fs::create_directory(p_cr_plus_ini_name);
		fs::path copy_ini = p_cr_plus_ini_name;
		copy_ini.append(ini_filename + ".inidem");
		fs::copy(ini_input_file, copy_ini.string(), fs::copy_options::overwrite_existing);
		std::cout << "copy [" << ini_input_file << "]" << std::endl;

		fs::path p_output_1 = p_cr_plus_ini_name;
		p_output_1.append("OUTPUT");
		if (fs::exists(p_output_1)) fs::remove_all(p_output_1);
		fs::create_directory(p_output_1);
		std::cout << "create [OUTPUT] directory" << std::endl;

		fs::path p_ptcl_stay = p_output_1;
		p_ptcl_stay.append(ini_filename + "_stay_ptcl.csv");
		ptcl_stay_last_file = p_ptcl_stay.string();

		fs::path p_ptcl_out = p_output_1;
		p_ptcl_out.append(ini_filename + "_out_ptcl.csv");
		ptcl_out_last_file = p_ptcl_out.string();

		fs::path p_time = p_output_1;
		p_time.append(ini_filename + "_time.txt");
		time_file = p_time.string();

		if (is_ptcl_output)
		{
			fs::path p_output_2 = p_cr_plus_ini_name;
			p_output_2.append("OUTPUT_PTCL");
			if (fs::exists(p_output_2)) fs::remove_all(p_output_2);
			fs::create_directory(p_output_2);
			std::cout << "create [OUTPUT_PTCL] directory" << std::endl;

			fs::path p_ptcl_output = p_output_2;
			p_ptcl_output.append("_");
			ptcl_file = p_ptcl_output.string();
			scope.ShowVariable(ptcl_file, "ptcl_interval_file");
		}

		if (is_avs_output)
		{
			fs::path p_avs = p_output_1;
			p_avs.append(ini_filename + "_ptcl.mgf");
			avs_file = p_avs.string();
			scope.ShowVariable(avs_file, "avs_file");
		}

		if (is_vtk_output)
		{
			fs::path p_output_3 = p_cr_plus_ini_name;
			p_output_3.append("OUTPUT_VTK");
			if (fs::exists(p_output_3)) fs::remove_all(p_output_3);
			fs::create_directory(p_output_3);
			std::cout << "create [OUTPUT_VTK] directory" << std::endl;

			fs::path p_vtk = p_output_3;
			p_vtk.append("_");
			vtk_file = p_vtk.string();
			scope.ShowVariable(vtk_file, "vtk_file");
		}

		scope.ShowVariable(time_file, "time_file");
		scope.ShowVariable(ptcl_out_last_file, "ptcl_out_last_file");
		scope.ShowVariable(ptcl_stay_last_file, "ptcl_stay_last_file");
		std::cout << "finish setting path to output files" << std::endl;
		std::cout << std::endl;
	}


	void DEM::LoadExtData() noexcept
	{
		// load electrostatic field data
		if (is_ele_field)
		{
			ele_field.LoadField(ele_input_file, ele_amp);
		}
		
		// load magnetic field data
		if (is_mag_field)
		{
			mag_field.LoadField(mag_input_file, mag_amp);
		}
		
		// load partile data
		ptcl.LoadPtcl(ptcl_input_file);
		
		// load object data
		obj.LoadObj(dobj_input_file);
	}


	void DEM::LoadParamater() noexcept
	{
		//open the configuration file
		std::cout << "loading ini [" << ini_input_file << "] file..." << std::endl;
		std::ifstream fin(ini_input_file);
		if (!fin)	data_io.ErrorInput(ini_input_file);
		std::string line;

		while (!fin.eof() && getline(fin, line))
		{	
			std::vector<std::string> key;
			data_io.SplitLine(line, key);
			if (key.empty()) continue;

			SetVector3D(key, "Calculation_Domain_XYZ", domain);
			SetArray(key, "Enable_Periodic_Boundary_XYZ", is_periodic);
			SetVector3D(key, "Cell_Size_Ratio_Collision_XYZ", coll_cell_ratio);
			SetVector3D(key, "Cell_Size_Ratio_Long_Range_XYZ", lr_cell_ratio);
			SetFlagVariable(key, "Enable_Cell_Z", is_cell_z, cell_domain_z);

			SetVariable(key, "Time_Step", time_step);
			SetVariable(key, "Simulation_Time", simulation_time);
			SetVariable(key, "Output_Time", output_time);
			SetFlagUndefinedArray(key, "Enable_Particle_Info_Output", is_ptcl_output, ptcl_key);
			SetFlagUndefinedArray(key, "Enable_Vtk_Output", is_vtk_output, vtk_key);
			SetFlagVariable(key, "Enable_Avs_Output", is_avs_output, avs_key);
			SetFlagVariable(key, "Enable_Openmp", is_openmp, openmp_size);

			SetVector3D(key, "Gravitational_Acceleration_XYZ", gravity);
			SetVector3D(key, "Air_Velocity_XYZ", air_velo);
			SetFlagVariable(key, "Enable_Air_Drag", is_airdrag, viscosity);
			SetFlagVariable(key, "Enable_Cunningham_Correction", is_cunningham, mean_free_path);
			SetFlagVariable(key, "Enable_Reynolds_Correction", is_reynolds, density_air);
			SetArray(key, "Permittivity_Particle_Air_Object", permittivity_ptcl_air_obj);
			SetArray(key, "Permeability_Particle_Air", permeability_ptcl_air);
			
			SetArray(key, "Coef_Restitution_Object", coef_rest_obj);
			SetArray(key, "Coef_Restitution_Particle", coef_rest_ptcl);
			SetArray(key, "Friction_Coef_Object_Particle", fri_coef_obj_ptcl);
			SetFlagVariable(key, "Enable_Rolling_Friction", is_roll_fri, coef_roll_fri);
			SetFlagVariable(key, "Enable_Adhesion_Force", is_adh, coef_adh);
			SetVariable(key, "Enable_Image_Force_Object", is_image_force);
			
			SetVariable(key, "Enable_Electrostatic_Field", is_ele_field);
			if (is_ele_field)
			{
				SetVariable(key, "Enable_Coulomb_External_Field", is_Coulomb_ext_field);
				SetVariable(key, "Enable_Coulomb_Between_Particle", is_Coulomb_ptcl);
				SetVariable(key, "Enable_Dielectrophoresis_Force", is_diele_force);
				SetVariable(key, "Enable_Dielectrophoresis_Torque", is_diele_torque);
				SetVariable(key, "Enable_Dielectrophoresis_Interaction", is_diele_dipole);
				SetVariable(key, "Electrostatic_Effect_Update_Time", ele_update_time);
				SetArray(key, "Electrostatic_Field_Start_Stop", ele_start_stop_time);
				SetUndefinedPairArray(key, "Electrostatic_Field_On_Off", 
					ele_on_off_time, ele_input_file.size());
				SetUndefinedArray(key, "Electrostatic_Field_Amplification", 
					ele_amp, ele_input_file.size());
			}

			SetVariable(key, "Enable_Magnetic_Field", is_mag_field);
			if (is_mag_field)
			{
				SetVariable(key, "Enable_Magnetic_Force", is_mag_force);
				SetVariable(key, "Enable_Magnetic_Torque", is_mag_torque);
				SetVariable(key, "Enable_Magnetic_Interaction", is_mag_dipole);
				SetVariable(key, "Magnetic_Effect_Update_Time", mag_update_time);
				SetArray(key, "Magnetic_Field_Start_Stop", mag_start_stop_time);
				SetUndefinedPairArray(key, "Magnetic_Field_On_Off", 
					mag_on_off_time, mag_input_file.size());
				SetUndefinedArray(key, "Magnetic_Field_Amplification", 
					mag_amp, mag_input_file.size());
			}
		}
		
		total_step = static_cast<int>(simulation_time / time_step + 0.001);
		
		const double io = output_time / time_step + 0.001;
		intvl_output = io > 1.0 ? static_cast<int>(io) : 1;
		
		const double ie = ele_update_time / time_step + 0.001;
		intvl_ele = ie > 1.0 ? static_cast<int>(ie) : 1;

		const double im = mag_update_time / time_step + 0.001;
		intvl_mag = im > 1.0 ? static_cast<int>(im) : 1;

		coef_ele = 1.0 / (4.0 * PI * PERMITTIVITY * permittivity_ptcl_air_obj[1]);
		coef_mag = (PERMEABILITY * permeability_ptcl_air[1]) / (4.0 * PI);

		coef_ele_moment = (4.0 * PI * PERMITTIVITY)
			* ((permittivity_ptcl_air_obj[0] - permittivity_ptcl_air_obj[1])
				/ (permittivity_ptcl_air_obj[0] + permittivity_ptcl_air_obj[1] * 2.0));

		coef_mag_moment = (4.0 * PI / PERMEABILITY)
			* ((permeability_ptcl_air[0] - permeability_ptcl_air[1])
				/ (permeability_ptcl_air[0] + permeability_ptcl_air[1] * 2.0));

		coef_image = 1.0 / (4.0 * PI * PERMITTIVITY)
			* (permittivity_ptcl_air_obj[2] - permittivity_ptcl_air_obj[1])
			/ (permittivity_ptcl_air_obj[2] + permittivity_ptcl_air_obj[1]);

		hard_sphere_info.coef_rest_n_obj = coef_rest_obj[0] + 1.0;
		hard_sphere_info.coef_rest_t_obj = coef_rest_obj[1] + 1.0;
		hard_sphere_info.coef_fri_obj = fri_coef_obj_ptcl[0];
		hard_sphere_info.coef_rest_n_ptcl = coef_rest_ptcl[0] + 1.0;
		hard_sphere_info.coef_rest_t_ptcl = coef_rest_ptcl[1] + 1.0;
		hard_sphere_info.coef_fri_ptcl = fri_coef_obj_ptcl[1];
		hard_sphere_info.coef_adh = coef_adh;
		hard_sphere_info.coef_roll_fri = coef_roll_fri;
		
		air_drag_info.coef_newton = 0.44 * PI * density_air * 0.5;
		air_drag_info.coef_allen = 10.0 * PI * density_air * 0.5;
		air_drag_info.coef_stokes = 6.0 * PI * viscosity;
		const double reynolds_coef = 2.0 * density_air / viscosity;
		air_drag_info.reynolds_coef_sq = reynolds_coef * reynolds_coef;
		air_drag_info.mean_free_path = mean_free_path;
		air_drag_info.air_velo = air_velo;

		ShowValue();
		std::cout << "finish loading ini [" << ini_input_file << "] file" << std::endl;
		fin.close();
	}


	void DEM::ShowValue() const noexcept
	{
		scope.ShowVariable(domain, "domain");
		scope.ShowArray(is_periodic, "is_periodic");
		scope.ShowVariable(coll_cell_ratio, "coll_cell_ratio");
		scope.ShowVariable(lr_cell_ratio, "lr_cell_ratio");
		scope.ShowVariable(is_cell_z, "is_cell_size_sorting_z");
		scope.ShowVariable(cell_domain_z, "cell_domain_z");

		scope.ShowVariable(time_step, "time_step");
		scope.ShowVariable(simulation_time, "simulation_time");
		scope.ShowVariable(output_time, "output_time");
		scope.ShowVariable(is_ptcl_output, "is_ptcl_output");
		scope.ShowArray(ptcl_key, "ptcl_key");
		scope.ShowVariable(is_vtk_output, "is_vtk_output");
		scope.ShowArray(vtk_key, "vtk_key");
		scope.ShowVariable(is_avs_output, "is_avs_output");
		scope.ShowArray(avs_key, "avs_key");
		scope.ShowVariable(is_openmp, "is_openmp");
		scope.ShowVariable(openmp_size, "openmp_size");

		scope.ShowVariable(gravity, "gravity");
		scope.ShowVariable(air_velo, "air_velo");
		scope.ShowVariable(is_airdrag, "is_airdrag");
		scope.ShowVariable(viscosity, "viscosity");
		scope.ShowVariable(is_cunningham, "is_cunningham");
		scope.ShowVariable(mean_free_path, "mean_free_path");
		scope.ShowVariable(is_reynolds, "is_reynolds");
		scope.ShowVariable(density_air, "density_air");
		scope.ShowArray(permittivity_ptcl_air_obj, "permittivity_ptcl_air_obj");
		scope.ShowArray(permeability_ptcl_air, "permeability_ptcl_air");

		scope.ShowArray(coef_rest_obj, "coef_rest_obj");
		scope.ShowArray(coef_rest_ptcl, "coef_rest_ptcl");
		scope.ShowArray(fri_coef_obj_ptcl, "fri_coef_obj_ptcl");
		scope.ShowVariable(is_roll_fri, "is_roll_fri");
		scope.ShowVariable(coef_roll_fri, "coef_roll_fri");
		scope.ShowVariable(is_adh, "is_adh");
		scope.ShowVariable(coef_adh, "coef_adh");
		scope.ShowVariable(is_image_force, "is_image_force");

		scope.ShowVariable(is_ele_field, "is_ele_field");
		scope.ShowVariable(is_Coulomb_ext_field, "is_Coulomb_ext_field");
		scope.ShowVariable(is_Coulomb_ptcl, "is_Coulomb_ptcl");
		scope.ShowVariable(is_diele_force, "is_diele_force");
		scope.ShowVariable(is_diele_torque, "is_diele_torque");
		scope.ShowVariable(is_diele_dipole, "is_diele_dipole");
		scope.ShowVariable(ele_update_time, "ele_update_time");
		scope.ShowArray(ele_start_stop_time, "ele_start_stop_time");
		scope.ShowMultipleArray(ele_on_off_time, "ele_on_off_time");
		scope.ShowArray(ele_amp, "ele_amp");
		
		scope.ShowVariable(is_mag_field, "is_mag_field");
		scope.ShowVariable(is_mag_force, "is_mag_force");
		scope.ShowVariable(is_mag_torque, "is_mag_torque");
		scope.ShowVariable(is_mag_dipole, "is_mag_dipole");
		scope.ShowVariable(mag_update_time, "mag_update_time");
		scope.ShowArray(mag_start_stop_time, "mag_start_stop_time");
		scope.ShowMultipleArray(mag_on_off_time, "mag_on_off_time");
		scope.ShowArray(mag_amp, "mag_amp");
	}

}