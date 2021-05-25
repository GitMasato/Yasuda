#pragma once
#include <iostream>
#include <vector>
#include <algorithm>

#include "Vector.h"

namespace dem
{

	struct NeighbPtclInfo
	{
		int id;
		Vec3 offset;

		NeighbPtclInfo() noexcept = default;

		NeighbPtclInfo(const int ID, const Vec3& Offset) noexcept
			: id(ID), offset(Offset) {}
	};


	struct SurrCellInfo 
	{
		int id;
		Vec3 offset;

		SurrCellInfo() noexcept : id(0), offset(0.0, 0.0, 0.0) {}
		
		SurrCellInfo(const int ID, const Vec3& Offset) noexcept
			: id(ID), offset(Offset) {}

		bool operator < (const SurrCellInfo& another) const noexcept
		{
			return id < another.id;
		};
	};
	  

	struct CellComponent 
	{
		int id;
		std::vector<int> contained_ptcl;
		std::vector<SurrCellInfo> surr_cell;

		CellComponent() noexcept : id(0) {}
		

		int GetContainedPtclSize() const noexcept
		{
			return static_cast<int>(contained_ptcl.size());
		}


		int GetSurrCellSize() const noexcept
		{
			return static_cast<int>(surr_cell.size());
		}
			

	};

}