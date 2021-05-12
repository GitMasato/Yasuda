#include "Output.h"

namespace dem
{
void Output::CreateTimeFile(const int Total_Time_Step,
                                   const std::array<double, 2> Energy,
                                   const std::string& File) const noexcept
{
  std::ofstream fout;
  fout.open(File, std::ios_base::out | std::ios_base::trunc);
  if (!fout) data_io.ErrorOutput(File);

  fout << "DEM START (" << clock.GetLocalDate() << ")" << std::endl;
  std::cout << "DEM START (" << clock.GetLocalDate() << ")" << std::endl;
  fout.close();

  WriteTimeFile(0, Total_Time_Step, Energy, File);
}

void Output::WriteTimeFile(const int Time_Step,
                                  const int Total_Time_Step,
                                  const std::array<double, 2> Energy,
                                  const std::string& File) const noexcept
{
  std::ofstream fout;
  fout.open(File, std::ios_base::out | std::ios_base::app);
  if (!fout) data_io.ErrorOutput(File);
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

  fout << " step (" << Time_Step << "/" << Total_Time_Step << ")"
       << "  erapsed_time (" << clock.GetElapsedTime() / 3600 << ":"
       << (clock.GetElapsedTime() % 3600) / 60 << ":"
       << (clock.GetElapsedTime() % 3600) % 60 << ")"
       << "  energy[" << ptcl_stay_size << "] (" << Energy[0] << ","
       << Energy[1] << ")" << std::endl;

  std::cout << " step (" << Time_Step << "/" << Total_Time_Step << ")"
            << "  erapsed_time (" << clock.GetElapsedTime() / 3600 << ":"
            << (clock.GetElapsedTime() % 3600) / 60 << ":"
            << (clock.GetElapsedTime() % 3600) % 60 << ")"
            << "  energy[" << ptcl_stay_size << "] (" << Energy[0] << ","
            << Energy[1] << ")" << std::endl;

  fout.close();
}

void Output::WriteTimeFileEnd(const std::string& File) const noexcept
{
  std::ofstream fout;
  fout.open(File, std::ios_base::out | std::ios_base::app);
  if (!fout) data_io.ErrorOutput(File);

  fout << "DEM End (" << clock.GetLocalDate() << ")" << std::endl;
  std::cout << "DEM End (" << clock.GetLocalDate() << ")" << std::endl;

  fout.close();
}

void Output::CreatePtclStayFile(const std::string& filename) const
  noexcept
{
  std::ofstream fout;
  fout.open(filename, std::ios_base::out | std::ios_base::trunc);
  if (!fout) data_io.ErrorOutput(filename);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();

  fout << ptcl->GetPtclStaySize() << ","
       << "number" << std::endl;
  fout << "x(m),y(m),z(m),u(m/s),v(m/s),w(m/s),"
       << "omega1(rad/s),omega2(rad/s),omega3(rad/s),"
       << "sp_gravity(kg/m3),diameter(m),charge(C)" << std::endl;

  for (const auto& v : p)
  {
    const Vec3 pos = v.GetPos();
    const Vec3 velo = v.GetVelo();
    const Vec3 ang_velo = v.GetAngVelo();
    fout << pos.x << "," << pos.y << "," << pos.z << ",";
    fout << velo.x << "," << velo.y << "," << velo.z << ",";
    fout << ang_velo.x << "," << ang_velo.y << "," << ang_velo.z << ",";
    fout << v.GetDensity() << ",";
    fout << v.GetDiameter() << ",";
    fout << v.GetChg() << std::endl;
  }

  fout.close();
}

void Output::CreatePtclOutFile(const std::string& filename) const
  noexcept
{
  std::ofstream fout;
  fout.open(filename, std::ios_base::out | std::ios_base::trunc);
  if (!fout) data_io.ErrorOutput(filename);

  ptcl->SortPtclOut();
  std::vector<ParticleInfo>& p = ptcl->GetPtclOut();

  fout << ptcl->GetPtclOutSize() << ","
       << "number" << std::endl;
  fout << "x(m),y(m),z(m),u(m/s),v(m/s),w(m/s),"
       << "omega1(rad/s),omega2(rad/s),omega3(rad/s),"
       << "sp_gravity(kg/m3),diameter(m),charge(C)" << std::endl;

  for (const auto& v : p)
  {
    const Vec3 pos = v.GetPos();
    const Vec3 velo = v.GetVelo();
    const Vec3 ang_velo = v.GetAngVelo();
    fout << pos.x << "," << pos.y << "," << pos.z << ",";
    fout << velo.x << "," << velo.y << "," << velo.z << ",";
    fout << ang_velo.x << "," << ang_velo.y << "," << ang_velo.z << ",";
    fout << v.GetDensity() << ",";
    fout << v.GetDiameter() << ",";
    fout << v.GetChg() << std::endl;
  }

  fout.close();
}

void Output::CreatePtclFile(const double Elapsed_Time,
                                   const std::string& File,
                                   std::vector<std::string> Keys) const noexcept
{
  std::stringstream ss;
  ss << File << std::setw(5) << std::setfill('0')
     << int((Elapsed_Time * 1000) + 0.5);
  ss << "_ms.dem";
  std::string file = ss.str();

  std::ofstream fout;
  fout.open(file,
            std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
  if (!fout) data_io.ErrorOutput(file);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_size = ptcl->GetPtclStaySize();
  fout.write((const char*)&ptcl_size, sizeof(int));

  for (const auto& i : Keys)
  {
    if (i == "id") WriteID("id", fout, p);
    if (i == "cp") WriteCollNumPtcl("cp", fout, p);
    if (i == "co") WriteCollNumObj("co", fout, p);

    if (i == "ra") WriteRadius("ra", fout, p);
    if (i == "ma") WriteMass("ma", fout, p);
    if (i == "uc") WriteUnitCharge("uc", fout, p);
    if (i == "ch") WriteCharge("ch", fout, p);
    if (i == "ee") WriteEleDipoleEnergy("ee", fout);
    if (i == "me") WriteMagDipoleEnergy("me", fout);

    if (i == "po") WritePos("po", fout, p);
    if (i == "ve") WriteVelo("ve", fout, p);
    if (i == "av") WriteAngVelo("av", fout, p);
    if (i == "an") WriteAng("an", fout, p);
    if (i == "fo") WriteForce("fo", fout, p);
    if (i == "to") WriteTorque("to", fout, p);

    if (i == "ef") WriteCoulombExtForce("ef", fout, p);
    if (i == "pf") WriteCoulombPtclForce("pf", fout, p);
    if (i == "df") WriteDieleForce("df", fout, p);
    if (i == "if") WriteImageForce("if", fout, p);
    if (i == "mf") WriteMagForce("mf", fout, p);
    if (i == "af") WriteAdhForce("af", fout, p);

    if (i == "dt") WriteDieleTorque("dt", fout, p);
    if (i == "mt") WriteMagTorque("mt", fout, p);
  }

  fout.close();
}

void Output::WriteID(const char Key[3], std::ofstream& Fout,
                            std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, id");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const int v = i.GetID();
    Fout.write((const char*)&v, sizeof(int));
  }
}

void Output::WriteCollNumPtcl(const char Key[3], std::ofstream& Fout,
                                     std::vector<ParticleInfo>& P) const
  noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, collision_num_ptcl");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const unsigned long long v = i.GetCollNumPtcl();
    Fout.write((const char*)&v, sizeof(unsigned long long));
  }
}

void Output::WriteCollNumObj(const char Key[3], std::ofstream& Fout,
                                    std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, collision_num_object");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const unsigned long long v = i.GetCollNumObj();
    Fout.write((const char*)&v, sizeof(unsigned long long));
  }
}

void Output::WriteRadius(const char Key[3], std::ofstream& Fout,
                                std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, radius");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const float v = static_cast<float>(i.GetRadius());
    Fout.write((const char*)&v, sizeof(float));
  }
}

void Output::WriteMass(const char Key[3], std::ofstream& Fout,
                              std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, mass");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const float v = static_cast<float>(i.GetMass());
    Fout.write((const char*)&v, sizeof(float));
  }
}

void Output::WriteUnitCharge(const char Key[3], std::ofstream& Fout,
                                    std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, unit_charge");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const float v = static_cast<float>(i.GetUnitChgMicroCGram());
    Fout.write((const char*)&v, sizeof(float));
  }
}

void Output::WriteCharge(const char Key[3], std::ofstream& Fout,
                                std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, charge");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const float v = static_cast<float>(i.GetChg());
    Fout.write((const char*)&v, sizeof(float));
  }
}

void Output::WriteEleDipoleEnergy(const char Key[3],
                                         std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, ele_dipole_energy");
  Fout.write((const char*)&Key, sizeof(char[2]));
  std::vector<FieldInfoAtParticle>& f = ele_field->GetFieldInfoAtPtcl();

  for (const auto& i : f)
  {
    const float v = static_cast<float>(i.GetMomentEnergy());
    Fout.write((const char*)&v, sizeof(float));
  }
}

void Output::WriteMagDipoleEnergy(const char Key[3],
                                         std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, mag_dipole_energy");
  Fout.write((const char*)&Key, sizeof(char[2]));
  std::vector<FieldInfoAtParticle>& f = mag_field->GetFieldInfoAtPtcl();

  for (const auto& i : f)
  {
    const float v = static_cast<float>(i.GetMomentEnergy());
    Fout.write((const char*)&v, sizeof(float));
  }
}

void Output::WritePos(const char Key[3], std::ofstream& Fout,
                             std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, pos");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetPos();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteVelo(const char Key[3], std::ofstream& Fout,
                              std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, velo");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetVelo();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteAngVelo(const char Key[3], std::ofstream& Fout,
                                 std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, ang_velo");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetAngVelo();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteAng(const char Key[3], std::ofstream& Fout,
                             std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, ang");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetAng();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteForce(const char Key[3], std::ofstream& Fout,
                               std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, force");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetExtForce();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteTorque(const char Key[3], std::ofstream& Fout,
                                std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, torque");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetExtTorque();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteCoulombExtForce(const char Key[3], std::ofstream& Fout,
                                         std::vector<ParticleInfo>& P) const
  noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, coulomb_ext_force");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetCoulombExtForce();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteCoulombPtclForce(const char Key[3],
                                          std::ofstream& Fout,
                                          std::vector<ParticleInfo>& P) const
  noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, coulomb_ptcl_force");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetCoulombPtclForce();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteImageForce(const char Key[3], std::ofstream& Fout,
                                    std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, image_force");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetImageForce();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteDieleForce(const char Key[3], std::ofstream& Fout,
                                    std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, diele_force");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetDieleForce();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteMagForce(const char Key[3], std::ofstream& Fout,
                                  std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, mag_force");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetMagForce();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteAdhForce(const char Key[3], std::ofstream& Fout,
                                  std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, adh_force");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetAdhForce();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteDieleTorque(const char Key[3], std::ofstream& Fout,
                                     std::vector<ParticleInfo>& P) const
  noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, diele_torque");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetDieleTorque();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

void Output::WriteMagTorque(const char Key[3], std::ofstream& Fout,
                                   std::vector<ParticleInfo>& P) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, mag_torque");
  Fout.write((const char*)&Key, sizeof(char[2]));

  for (const auto& i : P)
  {
    const Float3 v = i.GetMagTorque();
    Fout.write((const char*)&v, sizeof(Float3));
  }
}

}  // namespace dem
