#pragma once
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <sstream>
#include <vector>
#include <omp.h>

#include "Cell.h"
#include "Collision.h"
#include "DataIO.h"
#include "Field.h"
#include "LongRange.h"
#include "LongRangeParticleList.h"
#include "Object.h"
#include "OutputAVS.h"
#include "OutputVTK.h"
#include "OutputGeneral.h"
#include "Particle.h"
#include "Scope.h"
#include "Vector.h"

namespace dem
{

	class DEM {

	private:

		Cell cell, lr_cell;
		Collision coll;
		DataIO data_io;
		Field ele_field, mag_field;
		LongRange ele_lr, mag_lr;
		LongRangeParticleList lr_ptcl_list;
		Object obj;
		OutputAVS avs;
		OutputVTK vtk;
		OutputGeneral output;
		Particle ptcl;
		Scope scope;

		// file name
		std::vector<std::string> mag_input_file, ele_input_file, dobj_input_file;
		std::string ptcl_input_file, ini_input_file;
		std::string ptcl_stay_last_file, ptcl_out_last_file;
		std::string time_file, ptcl_file;
		std::string avs_file;
		std::string vtk_file;

		// input parameter
		Vec3 domain;										// calculation_domain_length xyz (m)Å@Å@
		std::array<bool, 3> is_periodic;					// to apply periodic boundary condition xyz, yes:1 no:0
		Vec3 coll_cell_ratio;								// cell size ratio to particle size, to sort particles spatially for collision calculation xyz
		Vec3 lr_cell_ratio;									// cell size ratio to particle size, to sort particles spatially for dipole interactions
		bool is_cell_z;										// to sort particles spatially using cell in z-direction, yes:1 no:0
		double cell_domain_z;								// domain to be applied cells for collision or dipole calculations z (m)

		double time_step;									// time step  (s)
		double simulation_time;								// simulation time  (s)
		double output_time;									// time step to output the results  (s)
		bool is_ptcl_output;								// to output particle information in each time interval
		std::vector<std::string> ptcl_key;					// particle output info keywords
		bool is_vtk_output;									// to output vtk data
		std::vector<std::string> vtk_key;					// vtk data keywords
		bool is_avs_output;									// to output MicroAVS file, no color, particle charge, dominant force  yes:1 no:0
		std::string avs_key;								// avs data keywords
		bool is_openmp;										// to apply paralell computing by openmp, yes:1 no:0
		int openmp_size;									// thread number for openmp

		Vec3 gravity;										// gravitational acceleration xyz	(m/s2)
		Vec3 air_velo;										// air velocity xyz (m/s)
		bool is_airdrag;									// to apply air drag, yes:1 no:0
		double viscosity;									// air viscosity (kg/s/m2)	
		bool is_cunningham;									// to apply cunningham_correction, yes:1 no:0
		double mean_free_path;								// mean free path									// pressureÅ@(Pa)	
		bool is_reynolds;									// to apply reynolds_number effect to air drag,  yes:1 no:0
		double density_air;									// air density  (kg/m3)
		std::array<double, 3> permittivity_ptcl_air_obj;	// permittivity ratio
		std::array<double, 2> permeability_ptcl_air;		// permittivity ratio

		std::array<double, 2> coef_rest_obj;				// coefficient of restitution
		std::array<double, 2> coef_rest_ptcl;				// coefficient of restitution
		std::array<double, 2> fri_coef_obj_ptcl;			// friction_coefficient
		bool is_roll_fri;									// to apply rolling friction, yes:1 no:0	
		double coef_roll_fri;								// coefficient of rolling friction
		bool is_adh;										// to apply adhesion force, yes:1 no:0	
		double coef_adh;									// coefficient of adhesion force
		bool is_image_force;								// to apply image force, yes:1 no:0

		bool is_ele_field;									// to apply electrostatic field, yes:1 no:0
		bool is_Coulomb_ext_field;							// to apply coulomb force from external field, yes:1 no:0
		bool is_Coulomb_ptcl;								// to apply coulomb force from other particles, yes:1 no:0
		bool is_diele_force;								// to apply dielectrophoresis_force, yes:1 no:0
		bool is_diele_torque;								// to apply dielectrophoresis_torque, yes:1 no:0
		bool is_diele_dipole;								// to apply dipole interaction between particles, yes:1 no:0
		double ele_update_time;								// time step to update electrostatic force  (s)
		std::array<double, 2> ele_start_stop_time;			// when to start and stop applying electrostatic field (s)
		std::vector<std::array<double, 2>> ele_on_off_time;	// duration to apply or remove one electrostatic field (s)
		std::vector<double> ele_amp;						// amplification rate to electrostatic field strength
		
		bool is_mag_field;									// to apply magnetic field, yes:1 no:0
		bool is_mag_force;									// to apply magnetic force, yes:1 no:0
		bool is_mag_torque;									// to apply magnetic torque, yes:1 no:0
		bool is_mag_dipole;									// to apply electrostatic dipole interactions between particle, yes:1 no:0
		double mag_update_time;								// time step to update magnetic force  (s)	
		std::array<double, 2> mag_start_stop_time;			// when to start and stop applying electrostatic field (s)
		std::vector<std::array<double, 2>> mag_on_off_time;	// duration to apply or remove one electrostatic field (s)
		std::vector<double> mag_amp;						// amplification rate to electrostatic field strength
		
		// temporal value, coefficient value
		int total_step;
		int intvl_output;
		int intvl_ele;
		int intvl_mag;
		double coef_ele;
		double coef_mag;
		double coef_ele_moment;
		double coef_mag_moment;
		double coef_image;
		AirDragInfo air_drag_info;
		HardSphereInfo hard_sphere_info;

		constexpr static double PI = 3.14159;
		constexpr static double PERMITTIVITY = 8.85418e-12;
		constexpr static double PERMEABILITY = 1.25663e-6;

	public:

		DEM() noexcept
			: domain(0.0, 0.0, 0.0), is_periodic{ false, false, false }
			, coll_cell_ratio(0.0, 0.0, 0.0), lr_cell_ratio(0.0, 0.0, 0.0)
			, is_cell_z(false), cell_domain_z(0.0)
			, time_step(0.0), simulation_time(0.0), output_time(0.0)
			, is_ptcl_output(false), is_vtk_output(false), is_avs_output(false)
			, is_openmp(false), openmp_size(0)
			, gravity(0.0, 0.0, 0.0), air_velo(0.0, 0.0, 0.0)
			, is_airdrag(false), viscosity(0.0)
			, is_cunningham(false), mean_free_path(0.0)
			, is_reynolds(false), density_air(0.0)
			, permittivity_ptcl_air_obj{ 0.0, 0.0, 0.0 }
			, permeability_ptcl_air{ 0.0, 0.0 }
			, coef_rest_obj{ 0.0, 0.0 }, coef_rest_ptcl{ 0.0, 0.0 }
			, fri_coef_obj_ptcl{ 0.0, 0.0 }, is_roll_fri(false), coef_roll_fri(0.0)
			, is_adh(false), coef_adh(0.0), is_image_force(false)
			, is_ele_field(false), is_Coulomb_ext_field(false),is_Coulomb_ptcl(false)
			, is_diele_force(false), is_diele_torque(false), is_diele_dipole(false)							// to apply dipole interaction between particles, yes:1 no:0
			, ele_update_time(0.0), ele_start_stop_time{ 0.0, 0.0 }
			, is_mag_field(false), is_mag_force(false)
			, is_mag_torque(false), is_mag_dipole(false)
			, mag_update_time(0.0)	, mag_start_stop_time{ 0.0, 0.0 }
			, total_step(0), intvl_output(0)
			, intvl_ele(0), intvl_mag(0)
			, coef_ele(0.0), coef_mag(0.0)
			, coef_ele_moment(0.0), coef_mag_moment(0.0), coef_image(0.0)
			, air_drag_info(), hard_sphere_info() {}

		void RunDEM() noexcept;

		void CalcEleForce() noexcept;

		void CalcMagForce() noexcept;

		void CheckBoundary() noexcept;

		void SumForceTorque() noexcept;

		void AssociateObject() noexcept;

		void CreateCell() noexcept;

		void SetPathInput(const int& argc, char** argv) noexcept;

		void SetPathOutput() noexcept;

		void LoadExtData() noexcept;

		void LoadParamater() noexcept;

		void ShowValue() const noexcept;
			   

		template<typename T>
		void SetVariable(const std::vector<std::string>& Value, 
			const std::string& Key, T& Parameter) noexcept
		{
			if (Value[0] == Key)
			{
				if (Value.size() == 2)
				{
					data_io.TransformString(Value[1], Parameter);
				}
				else { data_io.ErrorInput(Key); }
			}
		}


		void SetVariable(const std::vector<std::string>& Value,
			const std::string& Key, std::string& Parameter) noexcept
		{
			if (Value[0] == Key)
			{
				if (Value.size() == 2)
				{
					data_io.TransformString(Value[1], Parameter);
				}
				else { data_io.ErrorInput(Key); }
			}
		}


		template<typename T>
		void SetFlagVariable(const std::vector<std::string>& Value, 
			const std::string& Key, bool& Frag, T& Parameter) noexcept
		{
			if (Value[0] == Key)
			{
				if (Value.size() == 3)
				{
					data_io.TransformString(Value[1], Frag);
					if (Frag) data_io.TransformString(Value[2], Parameter);
				}
				else { data_io.ErrorInput(Key); }
			}
		}


		void SetFlagVariable(const std::vector<std::string>& Value,
			const std::string& Key, bool& Frag, std::string& Parameter) noexcept
		{
			if (Value[0] == Key)
			{
				if (Value.size() == 3)
				{
					data_io.TransformString(Value[1], Frag);
					if (Frag) data_io.TransformString(Value[2], Parameter);
				}
				else { data_io.ErrorInput(Key); }
			}
		}


		template<typename T>
		void SetVector3D(const std::vector<std::string>& Value, 
			const std::string& Key, Vector3D<T>& Parameter) noexcept
		{
			if (Value[0] == Key)
			{
				if (Value.size() == 4)
				{
					data_io.TransformString(Value[1], Parameter.x);
					data_io.TransformString(Value[2], Parameter.y);
					data_io.TransformString(Value[3], Parameter.z);
				}
				else { data_io.ErrorInput(Key); }
			}
		}


		template<typename T>
		void SetArray(const std::vector<std::string>& Value, 
			const std::string& Key, T& Parameter) noexcept
		{
			if (Value[0] == Key)
			{
				const std::size_t values = Value.size();
				const std::size_t parameters = Parameter.size();

				if (values - 1 == parameters)
				{
					for (std::size_t i = 1; i < values; ++i)
					{
						data_io.TransformString(Value[i], Parameter[i - 1]);
					}
				}
				else { data_io.ErrorInput(Key); }
			}
		}


		template<typename T>
		void SetUndefinedArray(const std::vector<std::string>& Value, 
			const std::string& Key, T& Parameter, const std::size_t Check_Num) noexcept
		{
			if (Value[0] == Key)
			{
				const std::size_t values = Value.size();

				if ((values >= 2) && (values - 1 == Check_Num))
				{
					Parameter.resize(Check_Num);

					for (std::size_t i = 1; i < values; ++i)
					{
						data_io.TransformString(Value[i], Parameter[i - 1]);
					}
				}
				else { data_io.ErrorInput(Key); }
			}
		}

		
		template<typename T>
		void SetFlagUndefinedArray(const std::vector<std::string>& Value, 
			const std::string& Key, bool& Frag, T& Parameter) noexcept
		{
			if (Value[0] == Key)
			{
				const std::size_t values = Value.size();

				if (values >= 3)
				{
					data_io.TransformString(Value[1], Frag);

					if (Frag)
					{
						Parameter.resize(values - 2);

						for (std::size_t i = 2; i < values; ++i)
						{
							if (Frag) data_io.TransformString(Value[i], Parameter[i - 2]);
						}
					}
				}
				else { data_io.ErrorInput(Key); }
			}
		}


		template<typename T>
		void SetUndefinedPairArray(const std::vector<std::string>& Value, 
			const std::string& Key, T& Parameter, const std::size_t Check_Num) noexcept
		{
			if (Value[0] == Key)
			{
				const std::size_t values = Value.size();

				if ((values >= 3) && (values % 2 != 0) && (values / 2 == Check_Num))
				{
					Parameter.resize(Check_Num);

					for (std::size_t i = 0, num = Check_Num; i < num; ++i)
					{
						data_io.TransformString(Value[i * 2 + 1], Parameter[i][0]);
						data_io.TransformString(Value[i * 2 + 2], Parameter[i][1]);
					}
				}
				else { data_io.ErrorInput(Key); }
			}
		}


	};

}