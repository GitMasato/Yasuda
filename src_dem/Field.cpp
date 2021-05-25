#include "Field.h"

namespace dem
{

	void Field::ResizeArraySize() noexcept
	{
		field_at_ptcl.resize(ptcl->GetPtclStaySize());
	}


	void Field::ModifyPtcl() noexcept
	{
		if (ptcl->GetPtclStaySize() < GetFieldInfoAtPtclSize())
		{
			field_at_ptcl.clear();
			ResizeArraySize();
		}
	}

	   
	void Field::ShowValue() noexcept
	{
		for (const auto& n : field)
		{
			std::size_t i = &n - &field[0];
			std::cout << std::endl;
			std::cout << "field num : " << i << std::endl;

			const auto m = n.GetMeshSize();
			std::cout << "mesh_size (" << m[0] << "," << m[1] << "," << m[2] << ")" << std::endl;
			std::cout << "mesh_length (" << n.GetMeshLength() << ")" << std::endl;
			std::cout << "total_mesh_size (" << n.GetMeshSizeXYZ() << ")" << std::endl;
		}
	}
	

	void Field::LoadField(const std::vector<std::string>& File_Name,
		const std::vector<double>& Amp) noexcept
	{
		std::cout << "loading field data..." << std::endl;
		const int total_file_size = static_cast<int>(File_Name.size());
		field.resize(total_file_size);

		for (int i = 0; i < total_file_size; ++i)
		{
			std::cout << "open [" << File_Name[i] << "]" << std::endl;
			std::ifstream ifs(File_Name[i], std::ios_base::binary);
			if (!ifs) data_io.ErrorInput(File_Name[i]);

			double mesh_length;
			std::array<int, 3> mesh_size;
			ifs.read((char*)&mesh_length, sizeof(double));
			ifs.read((char*)&mesh_size, sizeof(std::array<int, 3>));

			field[i].SetMeshLength(mesh_length);
			field[i].SetMeshLengthInv(1.0 / mesh_length);
			field[i].SetMeshSize(mesh_size);

			if ((field[i].GetMeshSizeXY() <= 1)
				|| (field[i].GetMeshSizeXZ() <= 1)
				|| (field[i].GetMeshSizeYZ() <= 1))
			{
				data_io.ErrorInput("no field data");
			}

			std::vector<Vec3> f(field[i].GetMeshSizeXYZ());
			ifs.read((char*)f.data(), sizeof(Vec3) * field[i].GetMeshSizeXYZ());

			for (auto& v : f) v *= Amp[i];
			field[i].MoveSetField(f);
		}

		ShowValue();
		std::cout << "finish loading field data" << std::endl;
		std::cout << std::endl;

		/*
		// output template

		for (int i = 0; i < total_file_num; ++i)
		{
			ofstream ofs;
			ofs.open(Field_Filename[i], ios_base::out | ios_base::binary | ios_base::trunc);

			double mesh_length = 0.025;
			double mesh_length_inv = 1.0 / mesh_length;
			Vec3 domain(0.2, 0, 0.1);
			
			std::array<int, 3> mesh_size;
			mesh_size[0] = static_cast<int>(domain.x * mesh_length_inv + 1e-9) + 1;
			mesh_size[1] = static_cast<int>(domain.y * mesh_length_inv + 1e-9) + 1;
			mesh_size[2] = static_cast<int>(domain.z * mesh_length_inv + 1e-9) + 1;

			const int total_mesh_size = mesh_size[0] * mesh_size[0] * mesh_size[0];
			vector<Vec3> field(total_mesh_size, Vec3::UnitX());

			ofs.write((const char*)& mesh_length, sizeof(double));
			ofs.write((const char*)& mesh_size, sizeof(std::array<int, 3>));

			for (int j = 0; j < total_mesh_size; ++j)
			{
				ofs.write((const char*)& field[j], sizeof(Vec3));
			}
		}
		*/

	}


	void Field::SetTgtField(const double Time, const std::array<double, 2>& Start_Stop_Time,
		const std::vector<std::array<double, 2>>& On_Off_Time) noexcept
	{
		if ((Start_Stop_Time[0] <= Time) && (Time <= Start_Stop_Time[1]))
		{
			double one_cycle = 0.0;
			for (const auto& i : On_Off_Time) { one_cycle += (i[0] + i[1]); }
			const double time_from_start = fmod(Time - Start_Stop_Time[0], one_cycle);

			for (int i = 0, s = static_cast<int>(On_Off_Time.size()); i < s; ++i)
			{
				double start_field = 0.0;
				double stop_field = 0.0;

				for (int k = 0; k <= i; ++k)
				{
					start_field = On_Off_Time[k][0] + stop_field;
					stop_field += On_Off_Time[k][0] + On_Off_Time[k][1];
				}

				if (time_from_start <= start_field)
				{
					is_ext_field_on = true;
					tgt_field = i;
					return;
				}
				else if (time_from_start <= stop_field)
				{
					is_ext_field_on = false;
					tgt_field = -1;
					return;
				}
			}
		}
		else
		{
			is_ext_field_on = false;
			tgt_field = -1;
		}
	}


	void Field::SetPtclFieldInfo(const std::array<bool, 3>& is_Periodic) noexcept
	{
		std::array<int, 3> mesh_size = field[tgt_field].GetMeshSize();

		if (mesh_size[2] <= 1) SetPtclFieldInfo2DXY(is_Periodic);
		else if (mesh_size[1] <= 1)	SetPtclFieldInfo2DXZ(is_Periodic);
		else if (mesh_size[0] <= 1) SetPtclFieldInfo2DYZ(is_Periodic);
		else SetPtclFieldInfo3D(is_Periodic);
	}


	void Field::SetPtclFieldInfo2DXY(const std::array<bool, 3>& is_Periodic) noexcept
	{
		std::vector<ObjectInfo>& o = obj->GetFieldBase();
		const Vec2 o_post_pos = o[0].GetPostPosXY();
		
		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();
		const double m_double_inv = 0.5 * field[tgt_field].GetMeshLengthInv();

		if (is_Periodic[0] && is_Periodic[1])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXY() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);
				
				if (!(field[tgt_field].isPtclInField2DXY(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXY(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXYPeriodic(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXY(m_double_inv, left, right, front, back);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[0])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXY() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);
				
				if (!(field[tgt_field].isPtclInField2DXY(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXY(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXY(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXY(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXY(m_double_inv, left, right, front, back);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[1])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXY() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);
				
				if (!(field[tgt_field].isPtclInField2DXY(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXY(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXYPeriodic(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXY(m_double_inv, left, right, front, back);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXY() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField2DXY(int_pos)))
				{					
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXY(int_pos);				
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXY(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXY(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXY(m_double_inv, left, right, front, back);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}		
	}


	void Field::SetPtclFieldInfo2DXZ(const std::array<bool, 3>& is_Periodic) noexcept
	{
		std::vector<ObjectInfo>& o = obj->GetFieldBase();
		const Vec2 o_post_pos = o[0].GetPostPosXZ();

		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();
		const double m_double_inv = 0.5 * field[tgt_field].GetMeshLengthInv();
		
		if (is_Periodic[0] && is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField2DXZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXZ(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldXZPeriodic(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldXZPeriodic(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXZ(m_double_inv, left, right, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[0])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField2DXZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXZ(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldXZ(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldXZ(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXZ(m_double_inv, left, right, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);
				
				if (!(field[tgt_field].isPtclInField2DXZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXZ(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldXZPeriodic(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldXZPeriodic(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXZ(m_double_inv, left, right, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosXZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);
				
				if (!(field[tgt_field].isPtclInField2DXZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DXZ(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldXZ(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldXZ(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DXZ(m_double_inv, left, right, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}		
	}


	void Field::SetPtclFieldInfo2DYZ(const std::array<bool, 3>& is_Periodic) noexcept
	{
		std::vector<ObjectInfo>& o = obj->GetFieldBase();
		const Vec2 o_post_pos = o[0].GetPostPosYZ();

		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();
		const double m_double_inv = 0.5 * field[tgt_field].GetMeshLengthInv();
		
		if (is_Periodic[1] && is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosYZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField2DYZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DYZ(int_pos);
				const Vec3 front = field[tgt_field].GetYMinusFieldYZPeriodic(int_pos[0], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldYZPeriodic(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldYZPeriodic(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldYZPeriodic(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DYZ(m_double_inv, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[1])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosYZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField2DYZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DYZ(int_pos);
				const Vec3 front = field[tgt_field].GetYMinusFieldYZPeriodic(int_pos[0], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldYZPeriodic(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldYZ(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldYZ(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DYZ(m_double_inv, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosYZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField2DYZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DYZ(int_pos);
				const Vec3 front = field[tgt_field].GetYMinusFieldYZ(int_pos[0], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldYZ(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldYZPeriodic(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldYZPeriodic(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DYZ(m_double_inv, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec2 rel_post_pos = p[i].GetPostPosYZ() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField2DYZ(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh2DYZ(int_pos);
				const Vec3 front = field[tgt_field].GetYMinusFieldYZ(int_pos[0], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldYZ(int_pos[0], tgt);
				const Vec3 down = field[tgt_field].GetZMinusFieldYZ(int_pos[1], tgt);
				const Vec3 up = field[tgt_field].GetZPlusFieldYZ(int_pos[1], tgt);
				const auto fd = GetFieldDeriv2DYZ(m_double_inv, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}		
	}


	void Field::SetPtclFieldInfo3D(const std::array<bool, 3>& is_Periodic) noexcept
	{
		std::vector<ObjectInfo>& o = obj->GetFieldBase();
		const Vec3 o_post_pos = o[0].GetPostPos();

		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();
		const double m_double_inv = 0.5 * field[tgt_field].GetMeshLengthInv();
			
		if (is_Periodic[0] && is_Periodic[1] && is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3DPeriodic(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3DPeriodic(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[0] && is_Periodic[1])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3D(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3D(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[0] && is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXY(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXY(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3DPeriodic(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3DPeriodic(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[1] && is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3DPeriodic(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3DPeriodic(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[0])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusPeriodic(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusFieldPeriodic(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXY(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXY(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3D(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3D(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[1])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXYPeriodic(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3D(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3D(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else if (is_Periodic[2])
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXY(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXY(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3DPeriodic(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3DPeriodic(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}
		else
		{
			#pragma omp for schedule(dynamic)
			for (int i = 0; i < ptcl_stay_size; ++i)
			{
				const Vec3 rel_post_pos = p[i].GetPostPos() - o_post_pos;
				const auto int_pos = field[tgt_field].GetIntPos(rel_post_pos);

				if (!(field[tgt_field].isPtclInField3D(int_pos)))
				{
					field_at_ptcl[i].SetExtFieldZero();
					continue;
				}

				const int tgt = field[tgt_field].GetTgtMesh3D(int_pos);
				const Vec3 left = field[tgt_field].GetXMinusField(int_pos[0], tgt);
				const Vec3 right = field[tgt_field].GetXPlusField(int_pos[0], tgt);
				const Vec3 front = field[tgt_field].GetYMinusFieldXY(int_pos[1], tgt);
				const Vec3 back = field[tgt_field].GetYPlusFieldXY(int_pos[1], tgt);
				const Vec3 down = field[tgt_field].GetZMinusField3D(int_pos[2], tgt);
				const Vec3 up = field[tgt_field].GetZPlusField3D(int_pos[2], tgt);
				const auto fd = GetFieldDeriv3D(m_double_inv, left, right, front, back, down, up);
				field_at_ptcl[i].SetExtField(field[tgt_field].GetField(tgt), fd);
			}
		}		
	}


	void Field::SetMoment(const double Dipole_Coef) noexcept
	{
		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();
		
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const Vec3 f = field_at_ptcl[i].GetField();
			const Vec3 moment = (Dipole_Coef * p[i].GetRadiusCu()) * f;
			field_at_ptcl[i].SetMoment(moment);
		}
	}


	void Field::CalcCoulombExtForce() noexcept
	{
		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();

		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			p[i].SetCoulombExtForce(p[i].GetChg() * field_at_ptcl[i].GetField());
		}
	}
	   	 

	void Field::CalcDieleForce() noexcept
	{
		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();

		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const auto fd = field_at_ptcl[i].GetFieldDeriv();
			const Vec3 moment = field_at_ptcl[i].GetMoment();
			const Vec3 force(moment.Dot(fd[0]), moment.Dot(fd[1]), moment.Dot(fd[2]));
			p[i].SetDieleForce(force);
		}
	}


	void Field::CalcDieleTorque() noexcept
	{
		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();

		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const Vec3 f = field_at_ptcl[i].GetField();
			const Vec3 moment = field_at_ptcl[i].GetMoment();
			p[i].SetDieleTorque(moment.Cross(f));
		}
	}


	void Field::CalcMagForce() noexcept
	{
		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();

		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const auto fd = field_at_ptcl[i].GetFieldDeriv();
			const Vec3 moment = field_at_ptcl[i].GetMoment();
			const Vec3 force(moment.Dot(fd[0]), moment.Dot(fd[1]), moment.Dot(fd[2]));
			p[i].SetMagForce(force);
		}
	}


	void Field::CalcMagTorque() noexcept
	{
		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_stay_size = ptcl->GetPtclStaySize();

		#pragma omp for schedule(dynamic)
		for (int i = 0; i < ptcl_stay_size; ++i)
		{
			const Vec3 f = field_at_ptcl[i].GetField();
			const Vec3 moment = field_at_ptcl[i].GetMoment();
			p[i].SetMagTorque(moment.Cross(f));
		}
	}	

}
