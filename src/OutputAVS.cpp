#include "OutputAVS.h"

namespace dem
{
void OutputAVS::CreateAVS(const int Step, const Vec3& Gravity,
                          const Vec3& Domain, const std::string& Avs_File,
                          const std::string& Key) noexcept
{
  std::ofstream fout;
  fout.open(Avs_File, std::ios_base::out | std::ios_base::trunc);
  if (!fout) data_io.ErrorOutput(Avs_File);

  fout << "# Micro AVS Geom:2.00" << std::endl;
  fout << Step + 1 << std::endl;

  WriteHeader(0.0, fout);
  if (Key == "ch")
    WritePtclCharge(fout);
  else if (Key == "fo")
    WritePtclForce(fout, Gravity);
  else if (Key == "no")
    WritePtcl(fout);
  WriteObj(fout, Domain);
  fout.close();
}

void OutputAVS::WriteAVS(const double Time, const Vec3& Gravity,
                         const Vec3& Domain, const std::string& Avs_File,
                         const std::string& Key) noexcept
{
  std::ofstream fout;
  fout.open(Avs_File, std::ios_base::out | std::ios_base::app);
  if (!fout) data_io.ErrorOutput(Avs_File);

  WriteHeader(Time, fout);
  if (Key == "ch")
    WritePtclCharge(fout);
  else if (Key == "fo")
    WritePtclForce(fout, Gravity);
  else if (Key == "no")
    WritePtcl(fout);
  WriteObj(fout, Domain);
  fout.close();
}

void OutputAVS::WriteHeader(const double Time, std::ofstream& Fout) noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, avs_header");

  counter++;
  Fout << "step" << counter << " " << Time << "sec" << std::endl;
  Fout << "sphere" << std::endl;
  Fout << "DEM" << std::endl;
  Fout << "color" << std::endl;
  Fout << ptcl->GetPtclStaySize() << std::endl;
}

void OutputAVS::WritePtcl(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, ptcl_normal");
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  constexpr double color = 0.5;
  // gray_deep = 0.4, gray_light = 0.8, gray_high_light = 0.95, white = 1

  for (const auto& v : p)
  {
    const Vec3 p_pos = v.GetPos();

    Fout << p_pos.x << " " << p_pos.y << " " << p_pos.z << " " << v.GetRadius()
         << " " << color << " " << color << " " << color << " " << std::endl;
  }
}

void OutputAVS::WritePtclCharge(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, ptcl_charge");
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();

  for (const auto& v : p)
  {
    const Vec3 p_pos = v.GetPos();

    if (v.GetChg() >= 0)
    {
      Fout << p_pos.x << " " << p_pos.y << " " << p_pos.z << " "
           << v.GetRadius() << " 1 0 0 " << std::endl;
    }
    else
    {
      Fout << p_pos.x << " " << p_pos.y << " " << p_pos.z << " "
           << v.GetRadius() << " 0 0 1 " << std::endl;
    }
  }
}

void OutputAVS::WritePtclForce(std::ofstream& Fout, const Vec3& Gravity) const
  noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, ptcl_force");
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();

  for (const auto& v : p)
  {
    constexpr Vec3 adhesion(0.0, 1.0, 0.0);
    constexpr Vec3 coulomb_ptcl(1.0, 0.0, 1.0);
    constexpr Vec3 coulomb_ext(0.2, 0.2, 0.2);
    constexpr Vec3 diele(1.0, 0.5, 0.0);
    constexpr Vec3 gravity(0.0, 0.0, 0.0);
    constexpr Vec3 magnetic(0.1, 0.1, 1.0);
    constexpr Vec3 image(0.0, 1.0, 1.0);
    constexpr Vec3 no_force(0.5, 0.5, 0.5);

    constexpr std::array<Vec3, 8> color_vector{
      adhesion, coulomb_ptcl, coulomb_ext, diele,
      gravity,  magnetic,     image,       no_force};

    const double adhesion_f = v.GetAdhForce().LengthSq();
    const double coulomb_ptcl_f = v.GetCoulombPtclForce().LengthSq();
    const double coulomb_ext_f = v.GetCoulombExtForce().LengthSq();
    const double diele_f = v.GetDieleForce().LengthSq();
    const double gravity_f = (v.GetMass() * Gravity).LengthSq();
    const double magnetic_f = v.GetMagForce().LengthSq();
    const double image_f = v.GetImageForce().LengthSq();
    const double no_f = 1.0e-25;

    const std::array<double, 8> force_magnitude{
      adhesion_f, coulomb_ptcl_f, coulomb_ext_f, diele_f,
      gravity_f,  magnetic_f,     image_f,       no_f};

    const auto max =
      std::max_element(force_magnitude.begin(), force_magnitude.end());
    const int num =
      static_cast<int>(std::distance(force_magnitude.begin(), max));
    const Vec3 rgb = color_vector[num];
    const Vec3 p_pos = v.GetPos();

    Fout << p_pos.x << " " << p_pos.y << " " << p_pos.z << " " << v.GetRadius()
         << " " << rgb.x << " " << rgb.y << " " << rgb.z << " " << std::endl;
  }
}

void OutputAVS::WriteObj(std::ofstream& Fout, const Vec3& Domain) const noexcept
{
  if (obj->GetFloorBaseSize()) WriteFloor(Fout, Domain);
  if (obj->GetBoxSize()) WriteBox(Fout);
  if (obj->GetHollowBoxSize()) WriteHollowBox(Fout);
  if (obj->GetCylinderSize()) WriteCylinder(Fout);
  if (obj->GetHollowCylinderSize()) WriteHollowCylinder(Fout);
  if (obj->GetSphereSize()) WriteSphere(Fout);
  if (obj->GetHollowSphereSize()) WriteHollowSphere(Fout);
}

void OutputAVS::WriteFloor(std::ofstream& Fout, const Vec3& Domain) const
  noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, floor");

  std::vector<ObjectInfo>& o = obj->GetFloorBase();
  if (!o[0].isVisualOn()) return;

  const Vec3 o_pos = o[0].GetPos();

  Fout << "disjoint polygon" << std::endl;
  Fout << "floor" << std::endl;
  Fout << "facet" << std::endl;
  Fout << "color" << std::endl;
  Fout << 4 << std::endl;

  // gray_deep = 0.4, gray_light = 0.8, white = 1
  constexpr Vec3 rgb(0.8, 0.8, 0.8);
  constexpr double t = 0.5E-3;

  Fout << "4" << std::endl;
  Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z << " " << rgb.x << " "
       << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y << " " << o_pos.z << " " << rgb.x
       << " " << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y + Domain.y << " " << o_pos.z
       << " " << rgb.x << " " << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x << " " << o_pos.y + Domain.y << " " << o_pos.z << " " << rgb.x
       << " " << rgb.y << " " << rgb.z << " " << std::endl;

  Fout << "4" << std::endl;
  Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z << " " << rgb.x << " "
       << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z - t << " " << rgb.x << " "
       << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y << " " << o_pos.z - t << " "
       << rgb.x << " " << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y << " " << o_pos.z << " " << rgb.x
       << " " << rgb.y << " " << rgb.z << " " << std::endl;

  Fout << "4" << std::endl;
  Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z << " " << rgb.x << " "
       << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z - t << " " << rgb.x << " "
       << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x << " " << o_pos.y + Domain.y << " " << o_pos.z - t << " "
       << rgb.x << " " << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x << " " << o_pos.y + Domain.y << " " << o_pos.z << " " << rgb.x
       << " " << rgb.y << " " << rgb.z << " " << std::endl;

  Fout << "4" << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y << " " << o_pos.z << " " << rgb.x
       << " " << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y << " " << o_pos.z - t << " "
       << rgb.x << " " << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y + Domain.y << " " << o_pos.z - t
       << " " << rgb.x << " " << rgb.y << " " << rgb.z << " " << std::endl;
  Fout << o_pos.x + Domain.x << " " << o_pos.y + Domain.y << " " << o_pos.z
       << " " << rgb.x << " " << rgb.y << " " << rgb.z << " " << std::endl;
}

void OutputAVS::WriteBox(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, box");

  std::vector<BoxInfo>& o = obj->GetBox();

  for (const auto& v : o)
  {
    if (!v.isVisualOn()) continue;

    Fout << "disjoint polygon" << std::endl;
    Fout << "box" << std::endl;
    Fout << "facet" << std::endl;
    Fout << "color" << std::endl;
    Fout << 6 << std::endl;

    const Vec3 o_pos_1 = v.GetPos();
    const Vec3 o_pos_2 = o_pos_1 + v.GetDim();

    // gray_deep = 0.4, gray_light = 0.8, white = 1
    constexpr Vec3 rgb(0.8, 0.8, 0.8);

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb.x
         << " " << rgb.y << " " << rgb.z << " " << std::endl;
  }
}

void OutputAVS::WriteHollowBox(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, hollow_box");

  std::vector<HollowBoxInfo>& o = obj->GetHollowBox();

  for (const auto& v : o)
  {
    if (!v.isVisualOn()) continue;

    Fout << "disjoint polygon" << std::endl;
    Fout << "hollow_box" << std::endl;
    Fout << "facet" << std::endl;
    Fout << "color" << std::endl;
    Fout << 12 << std::endl;

    const Vec3 o_pos_1 = v.GetPos();
    const Vec3 o_pos_2 = o_pos_1 + v.GetDim();
    const Vec3 o_pos_1_in = o_pos_1 + v.GetThickness();
    const Vec3 o_pos_2_in = o_pos_2 - v.GetThickness();

    // gray_deep = 0.4, gray_light = 0.8, white = 1
    constexpr Vec3 rgb_1(0.8, 0.8, 0.8);
    constexpr Vec3 rgb_2(0.9, 0.9, 0.9);

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_1.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_1.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
    Fout << o_pos_1.x << " " << o_pos_2.y << " " << o_pos_2.z << " " << rgb_1.x
         << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_1_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_1_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_2_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_2_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_1_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_1_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_1_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_1_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_1_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_2_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_2_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_1_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_1_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_2_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_2_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_1_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_2_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_2_in.y << " " << o_pos_1_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_2_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_2_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;

    Fout << "4" << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_1_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_1_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_2_in.x << " " << o_pos_2_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
    Fout << o_pos_1_in.x << " " << o_pos_2_in.y << " " << o_pos_2_in.z << " "
         << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z << " " << std::endl;
  }
}

void OutputAVS::WriteCylinder(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, cylinder");

  std::vector<CylinderInfo>& o = obj->GetCylinder();

  for (const auto& v : o)
  {
    if (!v.isVisualOn()) continue;
    Fout << "column" << std::endl;
    Fout << "cylinder" << std::endl;
    Fout << "dvertex_and_color" << std::endl;
    Fout << "64" << std::endl;
    Fout << 1 << std::endl;

    const Vec3 o_pos_1 = v.GetPos();
    const Vec3 o_pos_2 = o_pos_1 + (v.GetLength() * v.GetDirUnitVec());

    // gray_deep = 0.4, gray_light = 0.8, white = 1
    constexpr Vec3 rgb(0.8, 0.8, 0.8);

    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " "
         << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " "
         << v.GetOuterRadius() << " " << rgb.x << " " << rgb.y << " " << rgb.z
         << " " << std::endl;
  }
}

void OutputAVS::WriteHollowCylinder(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, hollow_cylinder");

  std::vector<HollowCylinderInfo>& o = obj->GetHollowCylinder();

  for (const auto& v : o)
  {
    if (!v.isVisualOn()) continue;

    constexpr int division_num = 24;
    constexpr Vec3 rgb_1(0.8, 0.8, 0.8);
    constexpr Vec3 rgb_2(0.0, 0.0, 0.0);

    Fout << "disjoint polygon" << std::endl;
    Fout << "hollow_cylinder" << std::endl;
    Fout << "facet" << std::endl;
    Fout << "color" << std::endl;
    Fout << division_num * 4 << std::endl;

    const Vec3 o_dir_u_vec = v.GetDirUnitVec();
    const Vec3 o_another_edge_corr = v.GetLength() * o_dir_u_vec;
    const Vec3 o_pos_1 = v.GetPos();
    const double o_inner_R = v.GetInnerRadius();
    const double o_outer_R = v.GetOuterRadius();

    const double division_rad = 360 * (PI / 180) / division_num;
    const Vec3 T_1_rot(0.0, 1.0, 0.0);
    Vec3 T_1 = v.RotatedReturn(T_1_rot);

    for (int i = 1; i <= division_num; ++i)
    {
      const double rad_T_2 = i * division_rad;
      const Vec3 T_2_rot(0.0, cos(rad_T_2), sin(rad_T_2));
      const Vec3 T_2 = v.RotatedReturn(T_2_rot);

      const Vec3 pos_1 = o_pos_1 + T_1 * o_inner_R;
      const Vec3 pos_2 = o_pos_1 + T_1 * o_outer_R;
      const Vec3 pos_3 = o_pos_1 + T_2 * o_outer_R;
      const Vec3 pos_4 = o_pos_1 + T_2 * o_inner_R;
      const Vec3 pos_1_ano = pos_1 + o_another_edge_corr;
      const Vec3 pos_2_ano = pos_2 + o_another_edge_corr;
      const Vec3 pos_3_ano = pos_3 + o_another_edge_corr;
      const Vec3 pos_4_ano = pos_4 + o_another_edge_corr;

      Fout << "4" << std::endl;
      Fout << pos_1.x << " " << pos_1.y << " " << pos_1.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_2.x << " " << pos_2.y << " " << pos_2.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_3.x << " " << pos_3.y << " " << pos_3.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_4.x << " " << pos_4.y << " " << pos_4.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

      Fout << "4" << std::endl;
      Fout << pos_1_ano.x << " " << pos_1_ano.y << " " << pos_1_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_2_ano.x << " " << pos_2_ano.y << " " << pos_2_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_3_ano.x << " " << pos_3_ano.y << " " << pos_3_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_4_ano.x << " " << pos_4_ano.y << " " << pos_4_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

      Fout << "4" << std::endl;
      Fout << pos_1.x << " " << pos_1.y << " " << pos_1.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_4.x << " " << pos_4.y << " " << pos_4.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_4_ano.x << " " << pos_4_ano.y << " " << pos_4_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_1_ano.x << " " << pos_1_ano.y << " " << pos_1_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

      Fout << "4" << std::endl;
      Fout << pos_2.x << " " << pos_2.y << " " << pos_2.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_3.x << " " << pos_3.y << " " << pos_3.z << " " << rgb_1.x
           << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_3_ano.x << " " << pos_3_ano.y << " " << pos_3_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;
      Fout << pos_2_ano.x << " " << pos_2_ano.y << " " << pos_2_ano.z << " "
           << rgb_1.x << " " << rgb_1.y << " " << rgb_1.z << " " << std::endl;

      T_1 = T_2;
    }

    /*
    Fout << "circle" << std::endl;
    Fout << "hollow_cylinder_c" << std::endl;
    Fout << "empty" << std::endl;
    Fout << "normal_and_color" << std::endl;
    Fout << division_num << std::endl;
    Fout << 4 << std::endl;

    const Vec3 o_pos_2 = o_pos_1 + o_another_edge_corr;

    Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << " "
            << o_dir_u_vec.x << " " << o_dir_u_vec.y << " " << o_dir_u_vec.z <<
    " "
            << o_inner_R << " " << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z
    << std::endl; Fout << o_pos_1.x << " " << o_pos_1.y << " " << o_pos_1.z << "
    "
            << o_dir_u_vec.x << " " << o_dir_u_vec.y << " " << o_dir_u_vec.z <<
    " "
            << o_outer_R << " " << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z
    << std::endl;

    Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << " "
            << o_dir_u_vec.x << " " << o_dir_u_vec.y << " " << o_dir_u_vec.z <<
    " "
            << o_inner_R << " " << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z
    << std::endl; Fout << o_pos_2.x << " " << o_pos_2.y << " " << o_pos_2.z << "
    "
            << o_dir_u_vec.x << " " << o_dir_u_vec.y << " " << o_dir_u_vec.z <<
    " "
            << o_outer_R << " " << rgb_2.x << " " << rgb_2.y << " " << rgb_2.z
    << std::endl;
    */
  }
}

void OutputAVS::WriteSphere(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, sphere");

  std::vector<SphereInfo>& o = obj->GetSphere();

  for (const auto& v : o)
  {
    if (!v.isVisualOn()) continue;
    Fout << "sphere" << std::endl;
    Fout << "sphere" << std::endl;
    Fout << "color" << std::endl;
    Fout << 1 << std::endl;

    // gray_deep = 0.4, gray_light = 0.8, white = 1
    constexpr Vec3 rgb(0.8, 0.8, 0.8);
    const Vec3 o_pos = v.GetPos();

    Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z << " "
         << v.GetOuterRadius() << " " << rgb.x << " " << rgb.y << " " << rgb.z
         << " " << std::endl;
  }
}

void OutputAVS::WriteHollowSphere(std::ofstream& Fout) const noexcept
{
  if (!Fout) data_io.ErrorOutput("output avs, hollow_sphere");

  std::vector<HollowSphereInfo>& o = obj->GetHollowSphere();

  for (const auto& v : o)
  {
    if (!v.isVisualOn()) continue;
    Fout << "sphere" << std::endl;
    Fout << "hollow_sphere" << std::endl;
    Fout << "color" << std::endl;
    Fout << 2 << std::endl;

    // gray_deep = 0.4, gray_light = 0.8, white = 1
    constexpr Vec3 rgb(0.8, 0.8, 0.8);
    const Vec3 o_pos = v.GetPos();

    Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z << " "
         << v.GetInnerRadius() << " " << rgb.x << " " << rgb.y << " " << rgb.z
         << " " << std::endl;

    Fout << o_pos.x << " " << o_pos.y << " " << o_pos.z << " "
         << v.GetOuterRadius() << " " << rgb.x << " " << rgb.y << " " << rgb.z
         << " " << std::endl;
  }
}

}  // namespace dem
