#pragma once
#include <array>
#include <cmath>
#include <vector>   

#include "CellInfo.h"
#include "Particle.h"
#include "Object.h"
#include "Vector.h"
#include "Scope.h"

namespace dem
{
	class Cell
	{
	private:

		Particle* ptcl;
		Object* obj;
		DataIO data_io;
		Scope scope;
		std::vector<CellComponent> cell, link_info;
		std::vector<int> cell_ptcl_table;

		std::array<int, 3> cell_size;
		Vec3 cell_length;
		Vec3 cell_origin;
		constexpr static double PI = 3.14159;

	public:

		Cell() noexcept	
			: ptcl(nullptr), obj(nullptr)
			, cell_size({ 0, 0, 0 })
			, cell_length(0.0, 0.0, 0.0), cell_origin(0.0, 0.0, 0.0) {}
		

		int GetCellSizeXY() const noexcept
		{
			return cell_size[0] * cell_size[1];
		}


		int GetCellSizeXYZ() const noexcept
		{
			return cell_size[0] * cell_size[1] * cell_size[2];
		}


		int GetCellPtclTableSize() const noexcept 
		{
			return static_cast<int>(cell_ptcl_table.size());
		}
		

		Vec3 GetCellLength() const noexcept
		{
			return cell_length;
		}
		

		void AssociateObj(Object* Obj) noexcept 
		{
			obj = Obj; 
		}


		void AssociatePtcl(Particle* Ptcl) noexcept
		{
			ptcl = Ptcl;
		}


		void ShowValue() const noexcept;

		void ShowValueLinkInfo() const noexcept;

		void ClearContainedPtcl() noexcept;

		void ResizeCellPtclTable() noexcept;

		void AdaptToPtclSize() noexcept;
		
		bool isDuplicate(int ID, int Vs_ID) const noexcept;

		void CreateCell(
			const Vec3& Cell_Size_Ratio, const Vec3& Domain, 
			const bool is_Periodic_Z, const bool is_Cell_On_Z,
			const double Cell_Domain_Z) noexcept;

		void SetCellInfo() noexcept;

		void SetSurrCell(const int X, const int Y, const int Z) noexcept;

		void SetPeriodicEffect(
			const Vec3& Domain, const std::array<bool, 3>& is_Periodic) noexcept;

		void SetPeriodicEffect1DX(const double DomainX) noexcept;

		void SetPeriodicEffect1DY(const double DomainY) noexcept;

		void SetPeriodicEffect1DZ(const double DomainZ) noexcept;

		void SetPeriodicEffect2DXY(const double DomainX, const double DomainY) noexcept;

		void SetPeriodicEffect2DXZ(const double DomainX, const double DomainZ) noexcept;

		void SetPeriodicEffect2DYZ(const double DomainY, const double DomainZ) noexcept;

		void SetPeriodicEffect3D(
			const double DomainX, const double DomainY, const double DomainZ) noexcept;

		void SetLinkedSurrCell() noexcept;

		void SetPtclInCell(const bool is_Periodic_Z, const bool is_Cell_Z,
			const double Cell_Domain_Z) noexcept;

		void SetPtclInCellPeriodic3D() noexcept;

		void SetPtclInCell3D(const double Cell_Domain_Z) noexcept;

		void SetPtclInCell2D() noexcept;

		void GetNeighbPtclList(
			const int ID, std::vector<NeighbPtclInfo>& Ptcl_List) noexcept;


		int GetTgtMesh(const int X, const int Y, const int Z) const noexcept
		{
			if ((X >= 0) && (X < cell_size[0])
				&& (Y >= 0) && (Y < cell_size[1])
				&& (Z >= 0) && (Z < cell_size[2]))
			{
				return  GetCellSizeXY() * Z + cell_size[0] * Y + X;
			}
			else
			{
				return -1;
			}
		}

	};

}