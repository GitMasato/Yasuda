#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "Vector.h"
#include "DataIO.h"
#include "ObjectInfo.h"
#include "Scope.h"

namespace dem
{

	class Object {

	private:

		DataIO data_io;
		Scope scope;
		std::vector<ObjectInfo> floor_base;
		std::vector<ObjectInfo> field_base;
		std::vector<ObjectInfo> boundary_base;
		std::vector<BoxInfo> box;
		std::vector<HollowBoxInfo> hollow_box;
		std::vector<CylinderInfo> cylinder;
		std::vector<HollowCylinderInfo> hollow_cylinder;
		std::vector<SphereInfo> sphere;
		std::vector<HollowSphereInfo> hollow_sphere;
		constexpr static double PI = 3.14159;

	public:

		Object() noexcept = default;

		std::vector<ObjectInfo>& GetFloorBase() noexcept
		{
			return floor_base; 
		}

		std::vector<ObjectInfo>& GetFieldBase() noexcept
		{
			return field_base;
		}

		std::vector<ObjectInfo>& GetBoundaryBase() noexcept
		{
			return boundary_base;
		}

		std::vector<BoxInfo>& GetBox() noexcept
		{
			return box;
		}

		std::vector<HollowBoxInfo>& GetHollowBox() noexcept
		{
			return hollow_box; 
		}

		std::vector<CylinderInfo>& GetCylinder() noexcept
		{
			return cylinder;
		}

		std::vector<HollowCylinderInfo>& GetHollowCylinder() noexcept
		{
			return hollow_cylinder; 
		}

		std::vector<SphereInfo>& GetSphere() noexcept
		{ 
			return sphere;
		}

		std::vector<HollowSphereInfo>& GetHollowSphere() noexcept
		{
			return hollow_sphere; 
		}

		int GetFloorBaseSize() const noexcept 
		{ 
			return static_cast<int>(floor_base.size()); 
		}

		int GetFieldBaseSize() const noexcept 
		{ 
			return static_cast<int>(field_base.size()); 
		}

		int GetBoundaryBaseSize() const noexcept 
		{
			return static_cast<int>(boundary_base.size());
		}

		int GetBoxSize() const noexcept 
		{ 
			return static_cast<int>(box.size()); 
		}

		int GetHollowBoxSize() const noexcept
		{
			return static_cast<int>(hollow_box.size());
		}

		int GetCylinderSize() const noexcept
		{ 
			return static_cast<int>(cylinder.size());
		}

		int GetHollowCylinderSize() const noexcept 
		{ 
			return static_cast<int>(hollow_cylinder.size()); 
		}

		int GetSphereSize() const noexcept 
		{ 
			return static_cast<int>(sphere.size());
		}

		int GetHollowSphereSize() const noexcept
		{ 
			return static_cast<int>(hollow_sphere.size());
		}

		void ShowValue() noexcept;

		void ShowValueFloorBase() noexcept;

		void ShowValueFieldBase() noexcept;

		void ShowValueBoundaryBase() noexcept;

		void ShowValueBox() noexcept;

		void ShowValueHollowBox() noexcept;

		void ShowValueCylinder() noexcept;

		void ShowValueHollowCylinder() noexcept;

		void ShowValueSphere() noexcept;

		void ShowValueHollowSphere() noexcept;
			   
		void SetHStepVelo(const double Elapsed_Time) noexcept;

		void SetPostVelo(const double Elapsed_Time) noexcept;

		void SetPostPos(const double Elapsed_Time) noexcept;
		
		void LoadObj(const std::vector<std::string>& Obj_Filename) noexcept;

		void CreateObj(const std::string& Obj_Type,
			const std::vector<std::string>& Obj_Info) noexcept;

		void SetObjParameter(std::map<std::string, std::string>& Parameter, 
			const std::vector<std::string>& Obj_Info) noexcept;


		template<typename T>
		void SetBasicParameter(T& Obj, std::map<std::string, std::string>& P) noexcept
		{
			const Vec3 pos(
				data_io.StringToDouble(P["pos_x"]),
				data_io.StringToDouble(P["pos_y"]),
				data_io.StringToDouble(P["pos_z"]));
			Obj.SetPostPos(pos);
			Obj.SetStdPos(pos);
			Obj.SetPos(pos);

			Obj.SetVibAmp(data_io.StringToDouble(P["amplitude"]));
			Obj.SetVibFreq(data_io.StringToDouble(P["frequency"]));

			const Vec3 vec(
				data_io.StringToDouble(P["vib_x"]),
				data_io.StringToDouble(P["vib_y"]),
				data_io.StringToDouble(P["vib_z"]));
			Obj.SetVibUnitVec(vec);

			Obj.SetVibFlag(data_io.StringToBool(P["vibration_on"]));	
		}


		template<typename T>
		void SetOuterRadius(T& Obj, std::map<std::string, std::string>& P) noexcept
		{
			Obj.SetOuterRadius(data_io.StringToDouble(P["outer_radius"]));
		}


		template<typename T>
		void SetInnerRadius(T& Obj, std::map<std::string, std::string>& P) noexcept
		{
			Obj.SetInnerRadius(data_io.StringToDouble(P["inner_radius"]));
		}


		template<typename T>
		void SetCollFlagVisualFlag(T& Obj, std::map<std::string, std::string>& P) noexcept
		{
			Obj.SetCollFlag(data_io.StringToBool(P["collision_on"]));
			Obj.SetVisualFlag(data_io.StringToBool(P["visualization_on"]));
		}


		template<typename T>
		void SetBoxParameter(T& Obj, std::map<std::string, std::string>& P) noexcept
		{
			const Vec3 dim(
				data_io.StringToDouble(P["dim_x"]),
				data_io.StringToDouble(P["dim_y"]),
				data_io.StringToDouble(P["dim_z"]));
			Obj.SetDim(dim);
		}


		template<typename T>
		void SetCylinderParameter(T& Obj, std::map<std::string, std::string>& P) noexcept
		{
			Obj.SetLength(data_io.StringToDouble(P["length"]));

			const Vec3 dir(
				data_io.StringToDouble(P["dir_x"]),
				data_io.StringToDouble(P["dir_y"]),
				data_io.StringToDouble(P["dir_z"]));
			Obj.SetDirUnitVec(dir);

			const Vec3 unit_vector(1.0, 0.0, 0.0);
			const Vec3 cross_vector = unit_vector.Cross(dir);
			const Vec3 axis = cross_vector.LengthSq() != 0.0
				? cross_vector.Normalized() : unit_vector;
			Obj.SetRotAxis(axis);

			double cos_alpha = unit_vector.Dot(dir);
			if (cos_alpha > 1.0)	cos_alpha = 1.0;
			if (cos_alpha < -1.0)	cos_alpha = -1.0;
			const double rad = acos(cos_alpha);

			Obj.SetRotRad(rad);
			Obj.SetSinPlus(sin(rad));
			Obj.SetSinMinus(sin(-rad));
			Obj.SetCosPlus(cos(rad));
			Obj.SetCosMinus(cos(-rad));
		}


		template<typename T>
		void ShowValue(const T& Obj) noexcept
		{	
			const Vec3 pos = Obj.GetPos();
			const Vec3 velo = Obj.GetVelo();
			const Vec3 vib_u_vec = Obj.GetVibUnitVec();

			std::cout << "ID (" << Obj.GetID() << ")" << std::endl;
			std::cout << "Pos (" << pos.x << "," << pos.y << "," << pos.z << ")" << std::endl;
			std::cout << "Velo (" << velo.x << "," << velo.y << "," << velo.z << ")" << std::endl;
			std::cout << "Coll_On (" << Obj.isCollOn() << ")" << std::endl;
			std::cout << "Vib_On (" << Obj.isVibOn() << ")" << std::endl;
			std::cout << "Vib_Amp (" << Obj.GetVibAmp() << ")" << std::endl;
			std::cout << "Vib_Freq (" << Obj.GetVibFreq() << ")" << std::endl;
			std::cout << "Vib_Unit_Vec (" << vib_u_vec.x << "," << vib_u_vec.y << "," << vib_u_vec.z << ")" << std::endl;
		}


	};

}