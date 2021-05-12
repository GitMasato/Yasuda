#pragma once
#include <algorithm>
#include <array>

#include "Vector.h"

namespace dem
{

	class ParticleInfo
	{
	private:

		int id;
		double radius, mass, mass_inv;
		double inertial_moment, inertial_moment_inv, chg;
		unsigned long long coll_num_ptcl, coll_num_obj;
		Vec3 post_pos, pos;
		Vec3 h_step_velo, velo;
		Vec3 ang, h_step_ang_velo, ang_velo;
		Vec3 post_ext_force, ext_force;
		Vec3 post_ext_torque, ext_torque;
		Vec3 coulomb_ext_force, diele_force, coulomb_ptcl_force, image_force;
		Vec3 mag_force, adh_force;
		Vec3 diele_torque, mag_torque;
		constexpr static double PI = 3.14159;

	public:

		ParticleInfo() noexcept
			: id(0)
			, radius(0.0), mass(0.0), mass_inv(0.0)
			, inertial_moment(0.0), inertial_moment_inv(0.0), chg(0.0)
			, coll_num_ptcl(0), coll_num_obj(0)
			, post_pos(0.0, 0.0, 0.0), pos(0.0, 0.0, 0.0)
			, h_step_velo(0.0, 0.0, 0.0), velo(0.0, 0.0, 0.0)
			, ang(0.0, 0.0, 0.0), h_step_ang_velo(0.0, 0.0, 0.0), ang_velo(0.0, 0.0, 0.0)
			, post_ext_force(0.0, 0.0, 0.0), ext_force(0.0, 0.0, 0.0)
			, post_ext_torque(0.0, 0.0, 0.0), ext_torque(0.0, 0.0, 0.0)
			, coulomb_ext_force(0.0, 0.0, 0.0), diele_force(0.0, 0.0, 0.0)
			, coulomb_ptcl_force(0.0, 0.0, 0.0), image_force(0.0, 0.0, 0.0)
			, mag_force(0.0, 0.0, 0.0), adh_force(0.0, 0.0, 0.0)
			, diele_torque(0.0, 0.0, 0.0), mag_torque(0.0, 0.0, 0.0) {}

		int GetID() const noexcept 
		{
			return id; 
		}

		unsigned long long GetCollNumPtcl() const noexcept
		{
			return coll_num_ptcl; 
		}

		unsigned long long GetCollNumObj() const noexcept
		{
			return coll_num_obj; 
		}

		double GetRadius() const noexcept 
		{
			return radius;
		}
		
		double GetRadiusSq() const noexcept 
		{
			return radius * radius; 
		}

		double GetRadiusCu() const noexcept 
		{
			return radius * radius * radius;
		}

		double GetDiameter() const noexcept
		{
			return radius * 2.0;
		}
		
		double GetDensity() const noexcept
		{
			constexpr double c = 1.0 / (4.0 / 3.0 * PI);
			return mass * c / GetRadiusCu();
		}

		double GetMass() const noexcept
		{
			return mass; 
		}

		double GetMassInv() const noexcept
		{
			return mass_inv; 
		}

		double GetInertialMoment() const noexcept 
		{
			return inertial_moment; 
		}

		double GetInertialMomentInv() const noexcept 
		{
			return inertial_moment_inv;
		}
		
		double GetChg() const noexcept 
		{
			return chg; 
		}

		double GetChgSq() const noexcept 
		{
			return chg * chg; 
		}

		double GetUnitChgMicroCGram() const noexcept
		{
			return chg + mass_inv * 1000.0;
		}
		
		Vec3 GetPostPos() const noexcept 
		{
			return post_pos;
		}

		Vec2 GetPostPosXY() const noexcept
		{
			return { post_pos.x, post_pos.y };
		}

		Vec2 GetPostPosXZ() const noexcept
		{
			return { post_pos.x, post_pos.z };
		}

		Vec2 GetPostPosYZ() const noexcept
		{
			return { post_pos.y, post_pos.z };
		}

		double GetPostPosX() const noexcept
		{
			return post_pos.x;
		}

		double GetPostPosY() const noexcept
		{
			return post_pos.y;
		}

		double GetPostPosZ() const noexcept
		{
			return post_pos.z;
		}

		Vec3 GetPos() const noexcept 
		{
			return pos; 
		}

		double GetPosZ() const noexcept
		{
			return pos.z;
		}

		Vec3 GetHStepVelo() const noexcept 
		{
			return h_step_velo;
		}

		double GetHStepVeloZ() const noexcept
		{
			return h_step_velo.z;
		}

		Vec3 GetVelo() const noexcept
		{
			return velo; 
		}

		Vec3 GetAng() const noexcept 
		{
			return ang;
		}

		Vec3 GetHStepAngVelo() const noexcept 
		{
			return h_step_ang_velo;
		}
		
		Vec3 GetAngVelo() const noexcept
		{
			return ang_velo;
		}

		Vec3 GetPostExtForce() const noexcept
		{
			return post_ext_force;
		}

		Vec3 GetExtForce() const noexcept 
		{
			return ext_force; 
		}

		Vec3 GetPostExtTorque() const noexcept 
		{
			return post_ext_torque;
		}

		Vec3 GetExtTorque() const noexcept 
		{
			return ext_torque; 
		}

		Vec3 GetCoulombExtForce() const noexcept 
		{
			return coulomb_ext_force; 
		}

		Vec3 GetDieleForce() const noexcept 
		{
			return diele_force; 
		}

		Vec3 GetCoulombPtclForce() const noexcept 
		{
			return coulomb_ptcl_force; 
		}

		Vec3 GetImageForce() const noexcept
		{
			return image_force; 
		}

		Vec3 GetMagForce() const noexcept 
		{
			return mag_force;
		}

		Vec3 GetAdhForce() const noexcept 
		{
			return adh_force;
		}
				
		Vec3 GetDieleTorque() const noexcept 
		{
			return diele_torque; 
		}

		Vec3 GetMagTorque() const noexcept 
		{
			return mag_torque; 
		}

		void SetID(const int v) noexcept
		{
			id = v;
		}

		void SetRadius(const double v) noexcept
		{
			radius = v;
		}
		
		void SetMass(const double v) noexcept
		{
			mass = v; 
		}

		void SetMassInv(const double v) noexcept
		{
			mass_inv = v;
		}

		void SetInertialMoment(const double v) noexcept
		{
			inertial_moment = v;
		}

		void SetInertialMomentInv(const double v) noexcept
		{
			inertial_moment_inv = v;
		}
		
		void SetChg(const double v) noexcept
		{
			chg = v;
		}
		
		void SetPostPos(const Vec3& v) noexcept
		{
			post_pos = v;
		}

		void SetPostPosX(const double v) noexcept
		{
			post_pos.x = v;
		}

		void SetPostPosY(const double v) noexcept
		{
			post_pos.y = v;
		}

		void SetPostPosZ(const double v) noexcept
		{
			post_pos.z = v;
		}

		void SetPos(const Vec3& v) noexcept 
		{
			pos = v; 
		}

		void SetHStepVelo(const Vec3& v) noexcept
		{
			h_step_velo = v;
		}

		void SetVelo(const Vec3& v) noexcept 
		{
			velo = v;
		}

		void SetAng(const Vec3& v) noexcept 
		{
			ang = v;
		}
				
		void SetHStepAngVelo(const Vec3& v) noexcept 
		{
			h_step_ang_velo = v;
		}
		
		void SetAngVelo(const Vec3& v) noexcept
		{
			ang_velo = v;
		}

		void SetPostExtForce(const Vec3& v) noexcept 
		{
			post_ext_force = v;
		}

		void SetExtForce(const Vec3& v) noexcept 
		{
			ext_force = v;
		}

		void SetPostExtTorque(const Vec3& v) noexcept 
		{
			post_ext_torque = v;
		}

		void SetExtTorque(const Vec3& v) noexcept 
		{
			ext_torque = v;
		}

		void SetCoulombExtForce(const Vec3& v) noexcept 
		{
			coulomb_ext_force = v;
		}

		void SetDieleForce(const Vec3& v) noexcept 
		{
			diele_force = v;
		}

		void SetCoulombPtclForce(const Vec3& v) noexcept 
		{
			coulomb_ptcl_force = v;
		}

		void SetImageForce(const Vec3& v) noexcept 
		{
			image_force = v;
		}

		void SetMagForce(const Vec3& v) noexcept
		{
			mag_force = v;
		}

		void SetAdhForce(const Vec3& v) noexcept 
		{
			adh_force = v;
		}
		
		void SetDieleTorque(const Vec3& v) noexcept
		{
			diele_torque = v;
		}

		void SetMagTorque(const Vec3& v) noexcept
		{
			mag_torque = v;
		}

		void AddCollNumPtcl() noexcept
		{
			coll_num_ptcl += 1;
		}

		void AddCollNumObj() noexcept
		{
			coll_num_obj += 1;
		}
		
		void AddPostPos(const Vec3& v) noexcept
		{
			post_pos += v;
		}
				
		void AddHStepVelo(const Vec3& v) noexcept
		{
			h_step_velo += v;
		}

		void AddAng(const Vec3& v) noexcept 
		{
			ang += v;
		}

		void AddHStepAngVelo(const Vec3& v) noexcept 
		{
			h_step_ang_velo += v;
		}
		
		void AddPostExtForce(const Vec3& v) noexcept
		{
			post_ext_force += v;
		}

		void AddPostExtTorque(const Vec3& v) noexcept
		{
			post_ext_torque += v;
		}

		void AddCoulombPtclForce(const Vec3& v) noexcept
		{
			#pragma omp atomic
			coulomb_ptcl_force.x += v.x;
			#pragma omp atomic
			coulomb_ptcl_force.y += v.y;
			#pragma omp atomic
			coulomb_ptcl_force.z += v.z;
		}

		void AddImageForce(const Vec3& v) noexcept
		{
			image_force += v;
		}
		
		void AddImageForceZ(const double v) noexcept
		{
			image_force.z += v;
		}
		
		void AddAdhForce(const Vec3& v) noexcept
		{
			#pragma omp atomic
			adh_force.x += v.x;
			#pragma omp atomic
			adh_force.y += v.y;
			#pragma omp atomic
			adh_force.z += v.z;
		}
		
		bool operator < (const ParticleInfo& another) const noexcept
		{
			return id < another.id;
		};

	};


	struct AirDragInfo
	{
		double coef_stokes;
		double coef_newton;
		double coef_allen;
		double reynolds_coef_sq;
		double mean_free_path;
		Vec3 air_velo;

		AirDragInfo() noexcept : coef_stokes(0.0), coef_newton(0.0), coef_allen(0.0),
			reynolds_coef_sq(0.0), mean_free_path(0.0), air_velo(0.0, 0.0, 0.0) {}
	};

}