#include "OutputVTK.h"

namespace dem
{

	void OutputVTK::CreateFile(const double Elapsed_time,
		const std::string& File_Name, std::vector<std::string> Keys) const noexcept
	{
		std::stringstream ss;
		ss << File_Name << std::setw(5) << std::setfill('0');
		ss << int((Elapsed_time * 1000) + 0.5) << "_ms";

		const std::string file = ss.str();
		ss << ".vtk";
		const std::string file_name = ss.str();

		std::ofstream fout;
		fout.open(file_name, std::ios_base::out | std::ios_base::trunc);
		if (!fout) data_io.ErrorOutput(file_name);

		std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
		const int ptcl_size = ptcl->GetPtclStaySize();

		fout << "# vtk DataFile Version 3.0" << std::endl;
		fout << file << std::endl;
		fout << "ASCII" << std::endl;
		fout << "DATASET UNSTRUCTURED_GRID" << std::endl << std::endl;

		WritePos(fout, p);
		WriteCellType(fout);
		fout << "POINT_DATA " << ptcl_size << std::endl;
		WriteRadius(fout, p);

		for (const auto& i : Keys)
		{
			if (i == "id") WriteID(fout, p);
			if (i == "cp") WriteCollNumPtcl(fout, p);
			if (i == "co") WriteCollNumObj(fout, p);

			if (i == "ma") WriteMass(fout, p);
			if (i == "uc") WriteUnitCharge(fout, p);
			if (i == "ch") WriteCharge(fout, p);
			if (i == "ee") WriteEleDipoleEnergy(fout);
			if (i == "me") WriteMagDipoleEnergy(fout);

			if (i == "ve") WriteVelo(fout, p);
			if (i == "av") WriteAngVelo(fout, p);
			if (i == "an") WriteAng(fout, p);
			if (i == "fo") WriteForce(fout, p);
			if (i == "to") WriteTorque(fout, p);

			if (i == "ef") WriteCoulombExtForce(fout, p);
			if (i == "pf") WriteCoulombPtclForce(fout, p);
			if (i == "df") WriteDieleForce(fout, p);
			if (i == "if") WriteImageForce(fout, p);
			if (i == "mf") WriteMagForce(fout, p);
			if (i == "af") WriteAdhForce(fout, p);

			if (i == "dt") WriteDieleTorque(fout, p);
			if (i == "mt") WriteMagTorque(fout, p);

		}

		fout.close();
	}


	void OutputVTK::WritePos(std::ofstream& Fout, std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, pos");
		const int ptcl_size = ptcl->GetPtclStaySize();
		Fout << "POINTS " << ptcl_size << " float" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetPos();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
		Fout << std::endl;
	}


	void OutputVTK::WriteCellType(std::ofstream& Fout) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, cell_type");
		const int ptcl_size = ptcl->GetPtclStaySize();
		Fout << "CELL_TYPES " << ptcl_size << std::endl;

		for (int i = 0; i < ptcl_size; ++i) { Fout << 1 << std::endl; }
		Fout << std::endl;
	}


	void OutputVTK::WriteID(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, id");
		Fout << "SCALARS id int" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P) { Fout << i.GetID() << std::endl; }
	}


	void OutputVTK::WriteCollNumPtcl(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, collision_num_ptcl");
		Fout << "SCALARS collision_num_ptcl int" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P) { Fout << i.GetCollNumPtcl() << std::endl; }
	}


	void OutputVTK::WriteCollNumObj(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, collision_num_object");
		Fout << "SCALARS collision_num_object int" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P) { Fout << i.GetCollNumObj() << std::endl; }
	}


	void OutputVTK::WriteRadius(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, radius");
		Fout << "SCALARS radius float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			Fout << static_cast<float>(i.GetRadius()) << std::endl;
		}
	}
	

	void OutputVTK::WriteMass(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, mass");
		Fout << "SCALARS mass float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			Fout << static_cast<float>(i.GetMass()) << std::endl;
		}
	}


	void OutputVTK::WriteUnitCharge(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, unit_charge");
		Fout << "SCALARS unit_charge float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			Fout << static_cast<float>(i.GetUnitChgMicroCGram()) << std::endl;
		}
	}


	void OutputVTK::WriteCharge(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, charge");
		Fout << "SCALARS charge float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			Fout << static_cast<float>(i.GetChg()) << std::endl;
		}
	}


	void OutputVTK::WriteEleDipoleEnergy(std::ofstream& Fout) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, ele_dipole_energy");
		Fout << "SCALARS ele_dipole_energy float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		std::vector<FieldInfoAtParticle>& f = ele_field->GetFieldInfoAtPtcl();

		for (const auto& i : f)
		{
			Fout << static_cast<float>(i.GetMomentEnergy()) << std::endl;
		}
	}


	void OutputVTK::WriteMagDipoleEnergy(std::ofstream& Fout) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, mag_dipole_energy");
		Fout << "SCALARS mag_dipole_energy float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		std::vector<FieldInfoAtParticle>& f = mag_field->GetFieldInfoAtPtcl();

		for (const auto& i : f)
		{
			Fout << static_cast<float>(i.GetMomentEnergy()) << std::endl;
		}
	}


	void OutputVTK::WriteVelo(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, velo");
		Fout << "VECTORS velo float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetVelo();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteAngVelo(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, ang_velo");
		Fout << "VECTORS ang_velo float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetAngVelo();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteAng(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, ang");
		Fout << "VECTORS ang float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetAng();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteForce(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, force");
		Fout << "VECTORS force float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetExtForce();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteTorque(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, torque");
		Fout << "VECTORS torque float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetExtTorque();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteCoulombExtForce(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, coulomb_ext_force");
		Fout << "VECTORS coulomb_ext_force float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetCoulombExtForce();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteCoulombPtclForce(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, coulomb_ptcl_force");
		Fout << "VECTORS coulomb_ptcl_force float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetCoulombPtclForce();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteImageForce(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, image_force");
		Fout << "VECTORS image_force float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetImageForce();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteDieleForce(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, diele_force");
		Fout << "VECTORS diele_force float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetDieleForce();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteMagForce(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, mag_force");
		Fout << "VECTORS mag_force float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetMagForce();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteAdhForce(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, adh_force");
		Fout << "VECTORS adh_force float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetAdhForce();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteDieleTorque(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, diele_torque");
		Fout << "VECTORS diele_torque float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetDieleTorque();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}


	void OutputVTK::WriteMagTorque(std::ofstream& Fout,
		std::vector<ParticleInfo>& P) const noexcept
	{
		if (!Fout) data_io.ErrorOutput("output vtk, mag_torque");
		Fout << "VECTORS mag_torque float" << std::endl;
		Fout << "LOOKUP_TABLE default" << std::endl;

		for (const auto& i : P)
		{
			const Float3 v = i.GetMagTorque();
			Fout << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}

}