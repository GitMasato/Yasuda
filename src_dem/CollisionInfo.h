#pragma once
#include <array>
#include <tuple>
#include <vector>

#include "Vector.h"

namespace dem
{

	struct CollisionInfo
	{
		// object: opponent object to collide withÅ@ 
		// 0: particle, 1: floor, 2: box, 3: hollow_box, 
		// 4: cylinder, 5: hollow_cylinder, 6: sphere, 6: hollow_sphere
		int obj_type, id_A, id_B;
		double coll_time, overlap;
		Vec3 vec_to_B, h_step_velo_B;

		CollisionInfo() noexcept = default;

		CollisionInfo(const int Type, const int ID_A, const int ID_B,
			const double Time, const double Overlap,
			const Vec3& Vec, const Vec3& H_Step_Velo_B) noexcept
			: obj_type(Type), id_A(ID_A), id_B(ID_B)
			, coll_time(Time), overlap(Overlap)
			, vec_to_B(Vec), h_step_velo_B(H_Step_Velo_B) {}


		bool operator < (const CollisionInfo& VS) const noexcept
		{
			constexpr double e = std::numeric_limits<double>::epsilon();

			return (std::abs(coll_time - VS.coll_time) > e) ? coll_time < VS.coll_time
				: (obj_type != VS.obj_type) ? obj_type < VS.obj_type
				: (id_A != VS.id_A) ? id_A < VS.id_A : id_B < VS.id_B;
		};


		std::tuple<int, double, Vec3, Vec3> GetCollObj() const noexcept
		{
			return { id_A, overlap, vec_to_B, h_step_velo_B };
		}


		std::tuple<int, int, double, Vec3> GetCollPtcl() const noexcept
		{
			return { id_A, id_B, overlap, vec_to_B };
		}

	};


	struct HardSphereInfo
	{
		double coef_rest_n_obj;
		double coef_rest_t_obj;
		double coef_fri_obj;
		double coef_rest_n_ptcl;
		double coef_rest_t_ptcl;
		double coef_fri_ptcl;
		double coef_adh;
		double coef_roll_fri;

		HardSphereInfo() noexcept :
			coef_rest_n_obj(0.0), coef_rest_t_obj(0.0), coef_fri_obj(0.0),
			coef_rest_n_ptcl(0.0), coef_rest_t_ptcl(0.0), coef_fri_ptcl(0.0),
			coef_adh(0.0), coef_roll_fri(0.0) {}

	};

}