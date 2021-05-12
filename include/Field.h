#pragma once
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "DataIO.h"
#include "FieldInfo.h"
#include "Object.h"
#include "Particle.h"
#include "Scope.h"
#include "Vector.h"

namespace dem
{
class Field
{
private:
  Particle* ptcl;
  Object* obj;
  DataIO data_io;
  Scope scope;

  bool is_ext_field_on;
  int tgt_field;
  std::vector<FieldInfo> field;
  std::vector<FieldInfoAtParticle> field_at_ptcl;

  constexpr static double PI = 3.14159;

public:
  Field() noexcept
    : ptcl(nullptr), obj(nullptr), is_ext_field_on(false), tgt_field(-1)
  {
  }

  void AssociatePtcl(Particle* Ptcl) noexcept { ptcl = Ptcl; }

  void AssociateObj(Object* Obj) noexcept { obj = Obj; }

  std::vector<FieldInfoAtParticle>& GetFieldInfoAtPtcl() noexcept
  {
    return field_at_ptcl;
  }

  bool isExtFieldOn() const noexcept { return is_ext_field_on; }

  int GetFieldSize() const noexcept { return static_cast<int>(field.size()); }

  int GetFieldInfoAtPtclSize() const noexcept
  {
    return static_cast<int>(field_at_ptcl.size());
  }

  void ResizeArraySize() noexcept;

  void ModifyPtcl() noexcept;

  void ShowValue() noexcept;

  void LoadField(const std::vector<std::string>& File_Name,
                 const std::vector<double>& Amp) noexcept;

  void SetTgtField(
    const double Time, const std::array<double, 2>& Start_Stop_Time,
    const std::vector<std::array<double, 2>>& On_Off_Time) noexcept;

  void SetPtclFieldInfo(const std::array<bool, 3>& is_Periodic) noexcept;

  void SetPtclFieldInfo2DXY(const std::array<bool, 3>& is_Periodic) noexcept;

  void SetPtclFieldInfo2DXZ(const std::array<bool, 3>& is_Periodic) noexcept;

  void SetPtclFieldInfo2DYZ(const std::array<bool, 3>& is_Periodic) noexcept;

  void SetPtclFieldInfo3D(const std::array<bool, 3>& is_Periodic) noexcept;

  void SetMoment(const double Dipole_Coef) noexcept;

  void CalcCoulombExtForce() noexcept;

  void CalcDieleForce() noexcept;

  void CalcDieleTorque() noexcept;

  void CalcMagForce() noexcept;

  void CalcMagTorque() noexcept;

  std::array<Vec3, 3> GetFieldDeriv2DXY(const double Mesh_Double_Inv,
                                        const Vec3& Left, const Vec3& Right,
                                        const Vec3& Front,
                                        const Vec3& Back) const noexcept
  {
    const Vec3 diff_x = (Right - Left) * Mesh_Double_Inv;
    const Vec3 diff_y = (Back - Front) * Mesh_Double_Inv;
    const Vec3 fdx(diff_x.x, diff_y.x, 0.0);
    const Vec3 fdy(diff_x.y, diff_y.y, 0.0);
    const Vec3 fdz(diff_x.z, diff_y.z, 0.0);
    return {fdx, fdy, fdz};
  }

  std::array<Vec3, 3> GetFieldDeriv2DXZ(const double Mesh_Double_Inv,
                                        const Vec3& Left, const Vec3& Right,
                                        const Vec3& Down, const Vec3& Up) const
    noexcept
  {
    const Vec3 diff_x = (Right - Left) * Mesh_Double_Inv;
    const Vec3 diff_z = (Up - Down) * Mesh_Double_Inv;
    const Vec3 fdx(diff_x.x, 0.0, diff_z.x);
    const Vec3 fdy(diff_x.y, 0.0, diff_z.y);
    const Vec3 fdz(diff_x.z, 0.0, diff_z.z);
    return {fdx, fdy, fdz};
  }

  std::array<Vec3, 3> GetFieldDeriv2DYZ(const double Mesh_Double_Inv,
                                        const Vec3& Front, const Vec3& Back,
                                        const Vec3& Down, const Vec3& Up) const
    noexcept
  {
    const Vec3 diff_y = (Back - Front) * Mesh_Double_Inv;
    const Vec3 diff_z = (Up - Down) * Mesh_Double_Inv;
    const Vec3 fdx(0.0, diff_y.x, diff_z.x);
    const Vec3 fdy(0.0, diff_y.y, diff_z.y);
    const Vec3 fdz(0.0, diff_y.z, diff_z.z);
    return {fdx, fdy, fdz};
  }

  std::array<Vec3, 3> GetFieldDeriv3D(const double Mesh_Double_Inv,
                                      const Vec3& Left, const Vec3& Right,
                                      const Vec3& Front, const Vec3& Back,
                                      const Vec3& Down, const Vec3& Up) const
    noexcept
  {
    const Vec3 diff_x = (Right - Left) * Mesh_Double_Inv;
    const Vec3 diff_y = (Back - Front) * Mesh_Double_Inv;
    const Vec3 diff_z = (Up - Down) * Mesh_Double_Inv;
    const Vec3 fdx(diff_x.x, diff_y.x, diff_z.x);
    const Vec3 fdy(diff_x.y, diff_y.y, diff_z.y);
    const Vec3 fdz(diff_x.z, diff_y.z, diff_z.z);
    return {fdx, fdy, fdz};
  }
};

}  // namespace dem
