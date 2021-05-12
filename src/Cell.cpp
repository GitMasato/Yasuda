#include "Cell.h"

namespace dem
{

	void Cell::ShowValue() const noexcept
	{
		scope.ShowArray(cell_size, "cell_size");
		scope.ShowVariable(cell_length, "cell_length");
		scope.ShowVariable(cell_origin, "cell_origin");

		for (const auto& c : cell)
		{
			for (const auto& r : c.surr_cell)
			{
				std::cout << "tgt, related cell ("
					<< c.id << "," << r.id << "," << r.offset
					<< ")" << std::endl;
			}
		}
	}


	void Cell::ShowValueLinkInfo() const noexcept
	{
		for (const auto& c : link_info)
		{
			for (const auto& r : c.surr_cell)
			{
				std::cout << "tgt, related cell ("
					<< c.id << "," << r.id << "," << r.offset
					<< ")" << std::endl;
			}
		}
	}


	void Cell::ClearContainedPtcl() noexcept
	{
		for (const auto& v : cell_ptcl_table) cell[v].contained_ptcl.clear();
	}


	void Cell::ResizeCellPtclTable() noexcept
	{
		cell_ptcl_table.clear();
		cell_ptcl_table.resize(ptcl->GetPtclStaySize(), 0);
	}


	void Cell::AdaptToPtclSize() noexcept
	{
		if (ptcl->GetPtclStaySize() < GetCellPtclTableSize())
		{
			for (const auto& v : cell_ptcl_table) cell[v].contained_ptcl.clear();			
			ResizeCellPtclTable();
		}
	}


	bool Cell::isDuplicate(int ID, int Vs_ID) const noexcept
	{
		if (Vs_ID == ID)
		{
			return true;
		}
		else
		{
			for(const auto& v : cell[ID].surr_cell)
			{
				if (Vs_ID == v.id) return true;
			}
			return false;
		}
	}


	void Cell::CreateCell(
		const Vec3& Cell_Size_Ratio, const Vec3& Domain,
		const bool is_Periodic_Z, const bool is_Cell_On_Z, 
		const double Cell_Domain_Z) noexcept
	{
		const double diameter_max = ptcl->GetMaxDiameter();
		cell_length = Cell_Size_Ratio * diameter_max;

		if (Cell_Size_Ratio.x < 1.0) cell_length.x = diameter_max;
		if (Cell_Size_Ratio.y < 1.0) cell_length.y = diameter_max;
		if (Cell_Size_Ratio.z < 1.0) cell_length.z = diameter_max;
		if (Domain.x / 3.0 < cell_length.x) cell_length.x = Domain.x / 3.0;
		if (Domain.y / 3.0 < cell_length.y) cell_length.y = Domain.y / 3.0;
		if (Domain.z / 3.0 < cell_length.z) cell_length.z = Domain.z / 3.0;		

		cell_size = (Domain / cell_length + 0.001).Int(); // 0.001: for rounding error
		for (auto v : cell_size) v += 1;
		cell_origin = ((cell_length * cell_size) - Domain) * 0.5;
		
		if (is_Periodic_Z)
		{
			if (!is_Cell_On_Z) data_io.ErrorInput("need to apply cell-sorting in z");
		}
		else if (is_Cell_On_Z)
		{			
			if (Cell_Size_Ratio.z < 1.0) cell_length.z = diameter_max;
			cell_size[2] = static_cast<int>(Cell_Domain_Z / cell_length.z + 0.001); // 0.001: for rounding error
			cell_origin.z = 0.0;
			
			if (Cell_Domain_Z < cell_length.z)
			{
				cell_length.z = Cell_Domain_Z;
				cell_size[2] = 1;
			}
		}
		else
		{			
			cell_size[2] = 1;	
			cell_origin.z = 0.0;
		}
			   		
		cell.resize(GetCellSizeXYZ());
		link_info.resize(GetCellSizeXYZ());
	}


	void Cell::SetCellInfo() noexcept
	{
		for (int k = 0; k < cell_size[2]; ++k)
		{
			for (int j = 0; j < cell_size[1]; ++j)
			{
				for (int i = 0; i < cell_size[0]; ++i)
				{
					const int id = GetTgtMesh(i, j, k);
					cell[id].id = id;
					link_info[id].id = id;
					SetSurrCell(i, j, k);
				}
			}
		}
	}


	void Cell::SetSurrCell(const int X, const int Y, const int Z) noexcept
	{
		const int tgt_id = GetTgtMesh(X, Y, Z);
		
		for (int k = 0; k < 3; ++k)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int i = 0; i < 3; ++i)
				{
					const int surr_id = GetTgtMesh(X - 1 + i, Y - 1 + j, Z - 1 + k);

					if (surr_id != -1 && surr_id != tgt_id)
					{
						cell[tgt_id].surr_cell.emplace_back(surr_id, Vec3::Zero());
					}
				}
			}
		}
	}


	void Cell::SetPeriodicEffect(
		const Vec3& Domain, const std::array<bool, 3>& is_Periodic) noexcept
	{
		if (is_Periodic[0] && is_Periodic[1] && is_Periodic[2])
		{
			std::cout << "periodic boundary 3D" << std::endl;
			SetPeriodicEffect3D(Domain.x, Domain.y, Domain.z);
		}
		else if (is_Periodic[0] && is_Periodic[1])
		{
			std::cout << "periodic boundary 2D XY" << std::endl;
			SetPeriodicEffect2DXY(Domain.x, Domain.y);
		}
		else if (is_Periodic[0] && is_Periodic[2])
		{
			std::cout << "periodic boundary 2D XZ" << std::endl;
			SetPeriodicEffect2DXZ(Domain.x, Domain.z);
		}
		else if (is_Periodic[1] && is_Periodic[2])
		{
			std::cout << "periodic boundary 2D YZ" << std::endl;
			SetPeriodicEffect2DYZ(Domain.y, Domain.z);
		}
		else if (is_Periodic[0])
		{
			std::cout << "periodic boundary 1D X" << std::endl;
			SetPeriodicEffect1DX(Domain.x);
		}
		else if (is_Periodic[1])
		{
			std::cout << "periodic boundary 1D Y" << std::endl;
			SetPeriodicEffect1DY(Domain.y);
		}
		else if (is_Periodic[2])
		{
			std::cout << "periodic boundary 1D Z" << std::endl;
			SetPeriodicEffect1DZ(Domain.z);
		}
	}


	void Cell::SetPeriodicEffect1DX(const double DomainX) noexcept
	{		
		const Vec3 id_0_offset(DomainX, 0.0, 0.0);
		const Vec3 id_1_offset(-DomainX, 0.0, 0.0);

		for (int k = 0; k < cell_size[2]; ++k)
		{
			for (int j = 0; j < cell_size[1]; ++j)
			{
				const int id_0 = GetTgtMesh(0, j, k);
				const int id_1 = GetTgtMesh(cell_size[0] - 1, j, k);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}
	}


	void Cell::SetPeriodicEffect1DY(const double DomainY) noexcept
	{		
		const Vec3 id_0_offset(0.0, DomainY, 0.0);
		const Vec3 id_1_offset(0.0, -DomainY, 0.0);

		for (int k = 0; k < cell_size[2]; ++k)
		{
			for (int i = 0; i < cell_size[0]; ++i)
			{
				const int id_0 = GetTgtMesh(i, 0, k);
				const int id_1 = GetTgtMesh(i, cell_size[1] - 1, k);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}
	}


	void Cell::SetPeriodicEffect1DZ(const double DomainZ) noexcept
	{		
		const Vec3 id_0_offset(0.0, 0.0, DomainZ);
		const Vec3 id_1_offset(0.0, 0.0, -DomainZ);

		for (int j = 0; j < cell_size[1]; ++j)
		{
			for (int i = 0; i < cell_size[0]; ++i)
			{
				const int id_0 = GetTgtMesh(i, j, 0);
				const int id_1 = GetTgtMesh(i, j, cell_size[2] - 1);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}
	}


	void Cell::SetPeriodicEffect2DXY(
		const double DomainX, const double DomainY) noexcept
	{
		SetPeriodicEffect1DX(DomainX);
		SetPeriodicEffect1DY(DomainY);

		// opposing corner 1
		{			
			const Vec3 id_0_offset(DomainX, DomainY, 0.0);
			const Vec3 id_1_offset(-DomainX, -DomainY, 0.0);

			for (int k = 0; k < cell_size[2]; ++k)
			{
				const int id_0 = GetTgtMesh(0, 0, k);
				const int id_1 = GetTgtMesh(cell_size[0] - 1, cell_size[1] - 1, k);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}

		// opposing corner 2
		{			
			const Vec3 id_0_offset(-DomainX, DomainY, 0.0);
			const Vec3 id_1_offset(DomainX, -DomainY, 0.0);

			for (int k = 0; k < cell_size[2]; ++k)
			{
				const int id_0 = GetTgtMesh(cell_size[0] - 1, 0, k);
				const int id_1 = GetTgtMesh(0, cell_size[1] - 1, k);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}
	}

	void Cell::SetPeriodicEffect2DXZ(
		const double DomainX, const double DomainZ) noexcept
	{
		SetPeriodicEffect1DX(DomainX);
		SetPeriodicEffect1DZ(DomainZ);

		// opposing corner 1
		{			
			const Vec3 id_0_offset(DomainX, 0.0, DomainZ);
			const Vec3 id_1_offset(-DomainX, 0.0, -DomainZ);

			for (int j = 0; j < cell_size[1]; ++j)
			{
				const int id_0 = GetTgtMesh(0, j, 0);
				const int id_1 = GetTgtMesh(cell_size[0] - 1, j, cell_size[2] - 1);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}

		// opposing corner 2
		{
			const Vec3 id_0_offset(-DomainX, 0.0, DomainZ);
			const Vec3 id_1_offset(DomainX, 0.0, -DomainZ);

			for (int j = 0; j < cell_size[1]; ++j)
			{
				const int id_0 = GetTgtMesh(cell_size[0] - 1, j, 0);
				const int id_1 = GetTgtMesh(0, j, cell_size[2] - 1);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}
	}

	void Cell::SetPeriodicEffect2DYZ(
		const double DomainY, const double DomainZ) noexcept
	{
		SetPeriodicEffect1DY(DomainY);
		SetPeriodicEffect1DZ(DomainZ);

		// opposing corner 1
		{			
			const Vec3 id_0_offset(0.0, DomainY, DomainZ);
			const Vec3 id_1_offset(0.0, -DomainY, -DomainZ);

			for (int i = 0; i < cell_size[0]; ++i)
			{
				const int id_0 = GetTgtMesh(i, 0, 0);
				const int id_1 = GetTgtMesh(i, cell_size[1] - 1, cell_size[2] - 1);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}

		// opposing corner 2
		{			
			const Vec3 id_0_offset(0.0, -DomainY, DomainZ);
			const Vec3 id_1_offset(0.0, DomainY, -DomainZ);

			for (int i = 0; i < cell_size[0]; ++i)
			{
				const int id_0 = GetTgtMesh(i, cell_size[1] - 1, 0);
				const int id_1 = GetTgtMesh(i, 0, cell_size[2] - 1);
				link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
				link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
			}
		}
	}


	void Cell::SetPeriodicEffect3D(
		const double DomainX, const double DomainY, const double DomainZ) noexcept
	{
		SetPeriodicEffect2DXY(DomainX, DomainY);
		SetPeriodicEffect2DXZ(DomainX, DomainZ);
		SetPeriodicEffect2DYZ(DomainY, DomainZ);

		// opposing corner 1
		{			
			const Vec3 id_0_offset(DomainX, DomainY, DomainZ);
			const Vec3 id_1_offset(-DomainX, -DomainY, -DomainZ);
			const int id_0 = GetTgtMesh(0, 0, 0);
			const int id_1 = GetTgtMesh(cell_size[0] - 1, cell_size[1] - 1, cell_size[2] - 1);
			link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
			link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
		}

		// opposing corner 2
		{
			const Vec3 id_0_offset(DomainX, -DomainY, DomainZ);
			const Vec3 id_1_offset(-DomainX, DomainY, -DomainZ);
			const int id_0 = GetTgtMesh(0, cell_size[1] - 1, 0);
			const int id_1 = GetTgtMesh(cell_size[0] - 1, 0, cell_size[2] - 1);
			link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
			link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
		}

		// opposing corner 3
		{
			const Vec3 id_0_offset(-DomainX, DomainY, DomainZ);
			const Vec3 id_1_offset(DomainX, -DomainY, -DomainZ);
			const int id_0 = GetTgtMesh(cell_size[0] - 1, 0, 0);
			const int id_1 = GetTgtMesh(0, cell_size[1] - 1, cell_size[2] - 1);
			link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
			link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
		}

		// opposing corner 4
		{
			const Vec3 id_0_offset(-DomainX, -DomainY, DomainZ);
			const Vec3 id_1_offset(DomainX, DomainY, -DomainZ);
			const int id_0 = GetTgtMesh(cell_size[0] - 1, cell_size[1] - 1, 0);
			const int id_1 = GetTgtMesh(0, 0, cell_size[2] - 1);
			link_info[id_0].surr_cell.emplace_back(id_1, id_1_offset);
			link_info[id_1].surr_cell.emplace_back(id_0, id_0_offset);
		}
	}


	void Cell::SetLinkedSurrCell() noexcept
	{
		std::vector<CellComponent> copy_cell;
		copy_cell = cell;
		
		// add surround cell of link cell of target cell
		for (int i = 0, c = GetCellSizeXYZ(); i < c; ++i)
		{
			for (int j = 0, lc = link_info[i].GetSurrCellSize(); j < lc; ++j)
			{
				const int link_cell_id = link_info[i].surr_cell[j].id;
				const int link_surr_cell_size = 
					copy_cell[link_cell_id].GetSurrCellSize();

				for (int k = 0; k < link_surr_cell_size; ++k)
				{
					const int link_surr_cell_id = 
						copy_cell[link_cell_id].surr_cell[k].id;
					if (isDuplicate(i, link_surr_cell_id)) continue;


					cell[i].surr_cell.emplace_back(
						link_surr_cell_id, link_info[i].surr_cell[j].offset);
				}
			}
		}


		copy_cell.clear();
		copy_cell = cell;

		// add linkcell of surround cell of target cell
		for (int i = 0, s = GetCellSizeXYZ(); i < s; ++i)
		{
			if (!(link_info[i].GetSurrCellSize())) continue;	

			for (int j = 0, sc = copy_cell[i].GetSurrCellSize(); j < sc; ++j)
			{
				const int surr_cell_id = copy_cell[i].surr_cell[j].id;
				const int surr_link_cell_size =
					link_info[surr_cell_id].GetSurrCellSize();

				for (int k = 0; k < surr_link_cell_size; ++k)
				{
					const int surr_link_cell_id = 
						link_info[surr_cell_id].surr_cell[k].id;

					if (isDuplicate(i, surr_link_cell_id)) continue;
					
					cell[i].surr_cell.emplace_back(
						surr_link_cell_id, 
						copy_cell[i].surr_cell[j].offset 
						+ link_info[surr_cell_id].surr_cell[k].offset);
				}
			}			
		}

		// add linkcell of target cell
		for (int i = 0, s = GetCellSizeXYZ(); i < s; ++i)
		{
			for (int j = 0, lc = link_info[i].GetSurrCellSize(); j < lc; ++j)
			{
				const int link_cell_id = link_info[i].surr_cell[j].id;

				if (isDuplicate(i, link_cell_id)) continue;

				cell[i].surr_cell.emplace_back(
					link_cell_id, 
					link_info[i].surr_cell[j].offset);
			}
		}


		copy_cell.clear();
		copy_cell = cell;

		// add linkcell of surround cell of target cell
		for (int i = 0, s = GetCellSizeXYZ(); i < s; ++i)
		{
			if (link_info[i].GetSurrCellSize()) continue;	

			for (int j = 0, sc = copy_cell[i].GetSurrCellSize(); j < sc; ++j)
			{
				const int surr_cell_id = copy_cell[i].surr_cell[j].id;
				const int surr_link_cell_size = 
					link_info[surr_cell_id].GetSurrCellSize();

				for (int k = 0; k < surr_link_cell_size; ++k)
				{
					const int surr_link_cell_id = 
						link_info[surr_cell_id].surr_cell[k].id;

					if (isDuplicate(i, surr_link_cell_id)) continue;

					cell[i].surr_cell.emplace_back(
						surr_link_cell_id,						
						link_info[surr_cell_id].surr_cell[k].offset);
				}
			}
		}

		for (int i = 0, s = GetCellSizeXYZ(); i < s; ++i)
		{
			std::sort(cell[i].surr_cell.begin(), cell[i].surr_cell.end());
		}

		link_info.clear();
		link_info.shrink_to_fit();
	}


	void Cell::SetPtclInCell(
		const bool is_Periodic_Z, const bool is_Cell_Z,
		const double Cell_Domain_Z) noexcept
	{		
		if (is_Periodic_Z)
		{
			if (is_Cell_Z)
			{
				SetPtclInCellPeriodic3D();
			}
			else 
			{
				data_io.ErrorInput("need to apply cell-sorting in z");
			}
		}
		else if (is_Cell_Z)
		{
			SetPtclInCell3D(Cell_Domain_Z);
		}
		else
		{
			SetPtclInCell2D();
		}
	}


	void Cell::SetPtclInCellPeriodic3D() noexcept
	{
		std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
		const Vec3 o_pos = o[0].GetPostPos() + cell_origin;
		const Vec3 cell_length_inv = 1.0 / cell_length;
		const int X = cell_size[0];
		const int XY = cell_size[0] * cell_size[1];

		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();
				
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const Vec3 rel_pos = p[i].GetPostPos() - o_pos;
			const auto tgt = (rel_pos * cell_length_inv).Int();
			const int id = XY * tgt[2] + X * tgt[1] + tgt[0];
			cell_ptcl_table[i] = id;

			#pragma omp critical
			cell[id].contained_ptcl.emplace_back(i);
		}	
	}


	void Cell::SetPtclInCell3D(const double Cell_Domain_Z) noexcept
	{
		std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
		const Vec3 o_pos = o[0].GetPostPos() + cell_origin;
		const Vec3 cell_length_inv = 1.0 / cell_length;
		const int X = cell_size[0];
		const int XY = cell_size[0] * cell_size[1];
		const int XY_Zminus = XY * (cell_size[2] - 1);

		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();

		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const Vec3 rel_pos = p[i].GetPostPos() - o_pos;
			const auto tgt = (rel_pos * cell_length_inv).Int();

			const int id = (rel_pos.z < 0.0) ? X * tgt[1] + tgt[0]
				: (Cell_Domain_Z < rel_pos.z) ? XY_Zminus + X * tgt[1] + tgt[0]
				: XY * tgt[2] + X * tgt[1] + tgt[0];

			cell_ptcl_table[i] = id;

			#pragma omp critical
			cell[id].contained_ptcl.emplace_back(i);
		}
	}


	void Cell::SetPtclInCell2D() noexcept
	{
		std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
		const Vec2 o_pos = o[0].GetPostPosXY() + cell_origin.XY();
		const Vec2 cell_length_inv = 1.0 / cell_length.XY();
		const int cell_x_size = cell_size[0];

		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();
		
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const Vec2 rel_pos = p[i].GetPostPosXY() - o_pos;
			const auto xy = (rel_pos * cell_length_inv).Int();
			const int id = cell_x_size * xy[1] + xy[0];
			cell_ptcl_table[i] = id;

			#pragma omp critical
			cell[id].contained_ptcl.emplace_back(i);
		}
	}


	void Cell::GetNeighbPtclList(
		const int ID, std::vector<NeighbPtclInfo>& Ptcl_List) noexcept
	{
		Ptcl_List.clear();
		const int tgt_cell_id = cell_ptcl_table[ID];

		for (const auto& v : cell[tgt_cell_id].contained_ptcl)
		{
			if (ID < v) Ptcl_List.emplace_back(v, Vec3::Zero());
		}

		for (const auto& c : cell[tgt_cell_id].surr_cell)
		{
			const int surr_cell_id = c.id;
			const Vec3 offset = c.offset;

			for (const auto& v : cell[surr_cell_id].contained_ptcl)
			{
				if (ID < v) Ptcl_List.emplace_back(v, offset);
			}
		}
	}
}