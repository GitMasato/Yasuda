#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "CollisionInfo.h"
#include "Particle.h"
#include "Object.h"
#include "Cell.h"
#include "Vector.h"
#include "Scope.h"

namespace dem
{

	class Collision
	{
	private:

		Particle* ptcl;
		Cell* cell;
		Object* obj;
		Scope scope;
		std::vector<CollisionInfo> coll_info;
		constexpr static double PI = 3.14159;

	public:

		Collision() noexcept : ptcl(nullptr), cell(nullptr), obj(nullptr) {}

		int GetCollInfoSize() const noexcept
		{
			return static_cast<int>(coll_info.size());
		}

		void AssociatePtcl(Particle* Ptcl) noexcept
		{
			ptcl = Ptcl;
		}

		void AssociateCell(Cell* Cell) noexcept
		{
			cell = Cell; 
		}

		void AssociateObj(Object* Obj) noexcept 
		{
			obj = Obj; 
		}

		void ClearCollOrder() noexcept;

		void DetectCollPtcl(const std::array<bool, 3>& is_Periodic, 
			const Vec3& Domain, const bool is_Openmp, const int Openmp_Size) noexcept;

		void DetectPtclColl(const int Reserve) noexcept;

		void DetectPtclCollPeriodicX(const Vec3& Domain, const int Reserve) noexcept;

		void DetectPtclCollPeriodicY(const Vec3& Domain, const int Reserve) noexcept;

		void DetectPtclCollPeriodicZ(const Vec3& Domain, const int Reserve) noexcept;

		void DetectPtclCollPeriodicXY(const Vec3& Domain, const int Reserve) noexcept;

		void DetectPtclCollPeriodicXZ(const Vec3& Domain, const int Reserve) noexcept;

		void DetectPtclCollPeriodicYZ(const Vec3& Domain, const int Reserve) noexcept;

		void DetectPtclCollPeriodicXYZ(const Vec3& Domain, const int Reserve) noexcept;
		
		void DetectCollObj(const bool is_Openmp, const int Openmp_Size) noexcept;

		void DetectCollObjImage(const bool is_Openmp, const int Openmp_Size, const double Image_Coef) noexcept;

		void DetectFloorColl(const int Reserve) noexcept;

		void DetectFloorCollImage(const double Image_Coef, const int Reserve) noexcept;

		void DetectBoxColl(const int Reserve) noexcept;

		void DetectBoxCollImage(const double Image_Coef, const int Reserve) noexcept;

		void DetectHollowBoxColl(const int Reserve) noexcept;

		void DetectHollowBoxCollImage(const double Image_Coef, const int Reserve) noexcept;

		void DetectCylinderColl(const int Reserve) noexcept;

		void DetectCylinderCollImage(const double Image_Coef, const int Reserve) noexcept;

		void DetectHollowCylinderColl(const int Reserve) noexcept;

		void DetectHollowCylinderCollImage(const double Image_Coef, const int Reserve) noexcept;

		void DetectSphereColl(const int Reserve) noexcept;

		void DetectSphereCollImage(const double Image_Coef, const int Reserve) noexcept;

		void DetectHollowSphereColl(const int Reserve) noexcept;

		void DetectHollowSphereCollImage(const double Image_Coef, const int Reserve) noexcept;

		void SortCollInfo() noexcept;

		void CalcColl(const HardSphereInfo& C) noexcept;

		void CalcCollAdh(const HardSphereInfo& C) noexcept;

		void CalcCollRoll(const HardSphereInfo& C) noexcept;

		void CalcCollAdhRoll(const HardSphereInfo& C) noexcept;


		int GetMaxPotentialCollNum() const noexcept
		{
			const Vec3 cell_len = cell->GetCellLength() / ptcl->GetMinDiameter();
			return static_cast<int>(cell_len.x * cell_len.y * cell_len.z);
		}


		void SetPeriodicEffect(double& Rel_Pos1D, const double Domein1D,
			const double H_Domain1D_Sq) const noexcept
		{
			if (Rel_Pos1D * Rel_Pos1D > H_Domain1D_Sq)
			{
				(Rel_Pos1D > 0.0) ? Rel_Pos1D -= Domein1D : Rel_Pos1D += Domein1D;
			}
		}


		std::pair<double, Vec3> GetOverlapVec(const double Dist_Sq,
			const double Min_dist, const Vec3& Rel_Pos) const noexcept
		{
			const double dist = sqrt(Dist_Sq);
			const double dist_inv = 1.0 / dist;
			return { Min_dist - dist, Rel_Pos * dist_inv };
		}
			   		 

		template<typename T>
		double GetCollTime(const T& Pre_Rel_Pos, const T& Pre_Rel_Velo,
			const double Dist_Sq_DIff) const noexcept
		{
			const double B = Pre_Rel_Pos.Dot(Pre_Rel_Velo);
			const double V = Pre_Rel_Velo.LengthSq();
			const double A = B * B - V * Dist_Sq_DIff;
			return (-B - sqrt(A)) / V;
		}		
		

		void ClampPos1D(bool& is_Inside_Box, double& O_Coll_Pt,
			const double Min, const double Max) const noexcept
		{
			if (O_Coll_Pt < Min)
			{
				O_Coll_Pt = Min;
				is_Inside_Box = false;
			}
			else if (O_Coll_Pt > Max)
			{
				O_Coll_Pt = Max;
				is_Inside_Box = false;
			}
		}


		bool ClampPos(Vec3& O_Coll_Pt, const Vec3& Min, const Vec3& Max) const noexcept
		{
			bool is_inside_box = true;
			ClampPos1D(is_inside_box, O_Coll_Pt.x, Min.x, Max.x);
			ClampPos1D(is_inside_box, O_Coll_Pt.y, Min.y, Max.y);
			ClampPos1D(is_inside_box, O_Coll_Pt.z, Min.z, Max.z);
			return is_inside_box;
		}


		std::pair<int, double> GetCloseDirDist(const Vec3& Ptcl,
			const Vec3& Box_Min, const Vec3& Box_Max) const noexcept
		{
			const std::array<double, 6> list
			{
				Ptcl.x - Box_Min.x, Box_Max.x - Ptcl.x,
				Ptcl.y - Box_Min.y, Box_Max.y - Ptcl.y,
				Ptcl.z - Box_Min.z, Box_Max.z - Ptcl.z
			};
			const auto min = std::min_element(list.begin(), list.end());
			const int dir = static_cast<int>(std::distance(list.begin(), min));
			return { dir, *min };
		}


		double GetDiffCloseDir(const int Dir, const Vec3& Ptcl,
			const Vec3& Box_Min, const Vec3& Box_Max) const noexcept
		{
			const std::array<double, 6> list
			{
				Ptcl.x - Box_Min.x, Box_Max.x - Ptcl.x,
				Ptcl.y - Box_Min.y, Box_Max.y - Ptcl.y,
				Ptcl.z - Box_Min.z, Box_Max.z - Ptcl.z,
			};
			return list[Dir];
		}


		constexpr Vec3 GetUNVecCloseDir(const int Dir) const noexcept
		{
			constexpr std::array<Vec3, 6> list
			{
				Vec3::UnitX(), -Vec3::UnitX(),
				Vec3::UnitY(), -Vec3::UnitY(),
				Vec3::UnitZ(), -Vec3::UnitZ()
			};
			return list[Dir];
		}
		

		double GetImageForceLen(const double Dist, const double Coef, const double ChgSq) const noexcept
		{
			const double image_dist = 2.0 * Dist;
			return Coef * ChgSq / (image_dist * image_dist);
		}
		

		std::pair<double, Vec3> GetJnLen(const Vec3& Pt_Velo, const Vec3& U_N_Vec, 
			const double Coef_Rest_N, const double M_Inv_Sum_Inv) const noexcept
		{
			const double pt_velo_n_len = Pt_Velo.Dot(U_N_Vec);
			const Vec3 pt_velo_t = Pt_Velo - (pt_velo_n_len * U_N_Vec);
			const double pt_velo_t_len_sq = pt_velo_t.LengthSq();
			const Vec3 u_t_vec = (pt_velo_t_len_sq <= 0.0) 
				? Vec3::Zero() : pt_velo_t / sqrt(pt_velo_t_len_sq);
			return { Coef_Rest_N * M_Inv_Sum_Inv * pt_velo_n_len, u_t_vec };
		}


		double GetJtLen(const Vec3& Pt_Velo, const Vec3& U_T_Vec, 
			const double Coef_Rest_T, const double M_Inv_Sum_Inv) const noexcept
		{
			constexpr double temp = 2.0 / 7.0;
			const double pt_velo_t_len = Pt_Velo.Dot(U_T_Vec);
			return temp * Coef_Rest_T * M_Inv_Sum_Inv * pt_velo_t_len;
		}
		

		Vec3 GetJ(const double Fri_Len, const double Jn_Len, const Vec3& U_N_Vec,
			const double Jt_Len, const Vec3& U_T_Vec) const noexcept
		{
			const Vec3 Jn = -Jn_Len * U_N_Vec;
			const Vec3 Jt = (Fri_Len > Jt_Len)
				? -Jt_Len * U_T_Vec : -Fri_Len * U_T_Vec;
			return Jn + Jt;
		}
		
		
		Vec3 GetRollFri(const Vec3& Ang_Velo, const double Jn_Len, const double Coef) const noexcept
		{	
			const double ang_velo_sq = Ang_Velo.LengthSq();
			if (ang_velo_sq <= 0.0) return Vec3::Zero();
			
			const double ang_velo_len = sqrt(ang_velo_sq);			
			const double roll_fri = std::abs(Jn_Len) * Coef * ang_velo_len;
			return (roll_fri >= ang_velo_len) ? Ang_Velo : roll_fri * Ang_Velo / ang_velo_len;
		}

	};

}