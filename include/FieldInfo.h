#pragma once
#include <array>
#include <vector>

#include "Vector.h"

namespace dem
{
class FieldInfoAtParticle
{
private:
  Vec3 moment, field_ptcl, field;
  std::array<Vec3, 3> field_deriv;

public:
  FieldInfoAtParticle() noexcept = default;

  constexpr void SetExtField(const Vec3& F,
                             const std::array<Vec3, 3>& Fd) noexcept
  {
    field = F;
    field_deriv = Fd;
  }

  constexpr void SetExtFieldZero() noexcept
  {
    field = Vec3::Zero();
    field_deriv[0] = Vec3::Zero();
    field_deriv[1] = Vec3::Zero();
    field_deriv[2] = Vec3::Zero();
  }

  void SetFieldAsTotalField() noexcept { field = GetFieldTotal(); }

  void SetMoment(const Vec3& F) noexcept { moment = F; }

  void SetFieldPtcl(const Vec3& F) noexcept { field_ptcl = F; }

  Vec3 GetField() const noexcept { return field; }

  Vec3 GetFieldTotal() const noexcept { return field + field_ptcl; }

  Vec3 GetMoment() const noexcept { return moment; }

  double GetMomentEnergy() const noexcept { return field.Dot(moment); }

  std::array<Vec3, 3> GetFieldDeriv() const noexcept { return field_deriv; }

  void AddFieldDeriv(const std::array<Vec3, 3>& Fd) noexcept
  {
    field_deriv[0] += Fd[0];
    field_deriv[1] += Fd[1];
    field_deriv[2] += Fd[2];
  }
};

class FieldInfo
{
private:
  std::array<int, 3> mesh_size;
  double mesh_len, mesh_len_inv;
  std::vector<Vec3> field;

public:
  FieldInfo() noexcept : mesh_size{0, 0, 0}, mesh_len(0.0), mesh_len_inv(0.0) {}

  std::array<int, 3> GetMeshSize() const noexcept { return mesh_size; }

  double GetMeshLength() const noexcept { return mesh_len; }

  double GetMeshLengthInv() const noexcept { return mesh_len_inv; }

  int GetMeshSizeXY() const noexcept { return mesh_size[0] * mesh_size[1]; }

  int GetMeshSizeXZ() const noexcept { return mesh_size[0] * mesh_size[2]; }

  int GetMeshSizeYZ() const noexcept { return mesh_size[1] * mesh_size[2]; }

  int GetMeshSizeXYZ() const noexcept
  {
    return mesh_size[0] * mesh_size[1] * mesh_size[2];
  }

  void SetMeshLength(const double v) noexcept { mesh_len = v; }

  void SetMeshLengthInv(const double v) noexcept { mesh_len_inv = v; }

  void SetMeshSize(const std::array<int, 3>& v) noexcept { mesh_size = v; }

  void MoveSetField(std::vector<Vec3>& v) noexcept { field = std::move(v); }

  std::array<int, 2> GetIntPos(const Vec2& Rel_Pos) const noexcept
  {
    const Vec2 pos_mesh = Rel_Pos * mesh_len_inv + 0.5;

    return {static_cast<int>(pos_mesh.x), static_cast<int>(pos_mesh.y)};
  }

  std::array<int, 3> GetIntPos(const Vec3& Rel_Pos) const noexcept
  {
    const Vec3 pos_mesh = Rel_Pos * mesh_len_inv + 0.5;

    return {static_cast<int>(pos_mesh.x), static_cast<int>(pos_mesh.y),
            static_cast<int>(pos_mesh.z)};
  }

  Vec3 GetField(const int Tgt) const noexcept { return field[Tgt]; }

  Vec3 GetXMinusField(const int X, const int Center) const noexcept
  {
    if (X == 0) return Vec3::Zero();

    const int n = Center - 1;
    return field[n];
  }

  Vec3 GetXMinusPeriodic(const int X, const int Center) const noexcept
  {
    const int n = (X == 0) ? Center + mesh_size[0] - 2 : Center - 1;
    return field[n];
  }

  Vec3 GetXPlusField(const int X, const int Center) const noexcept
  {
    if (X + 1 == mesh_size[0]) return Vec3::Zero();

    const int n = Center + 1;
    return field[n];
  }

  Vec3 GetXPlusFieldPeriodic(const int X, const int Center) const noexcept
  {
    const int n =
      (X + 1 == mesh_size[0]) ? Center - mesh_size[0] + 2 : Center + 1;
    return field[n];
  }

  Vec3 GetYMinusFieldXY(const int Y, const int Center) const noexcept
  {
    if (Y == 0) return Vec3::Zero();

    const int n = Center - mesh_size[0];
    return field[n];
  }

  Vec3 GetYMinusFieldXYPeriodic(const int Y, const int Center) const noexcept
  {
    const int n = (Y == 0) ? Center + GetMeshSizeXY() - (mesh_size[0] * 2)
                           : Center - mesh_size[0];
    return field[n];
  }

  Vec3 GetYPlusFieldXY(const int Y, const int Center) const noexcept
  {
    if (Y + 1 == mesh_size[1]) return Vec3::Zero();

    const int n = Center + mesh_size[0];
    return field[n];
  }

  Vec3 GetYPlusFieldXYPeriodic(const int Y, const int Center) const noexcept
  {
    const int n = (Y + 1 == mesh_size[1])
                    ? Center - GetMeshSizeXY() + (mesh_size[0] * 2)
                    : Center + mesh_size[0];
    return field[n];
  }

  Vec3 GetZMinusFieldXZ(const int Z, const int Center) const noexcept
  {
    if (Z == 0) return Vec3::Zero();

    const int n = Center - mesh_size[0];
    return field[n];
  }

  Vec3 GetZMinusFieldXZPeriodic(const int Z, const int Center) const noexcept
  {
    const int n = (Z == 0) ? Center + GetMeshSizeXZ() - (mesh_size[0] * 2)
                           : Center - mesh_size[0];
    return field[n];
  }

  Vec3 GetZPlusFieldXZ(const int Z, const int Center) const noexcept
  {
    if (Z + 1 == mesh_size[2]) return Vec3::Zero();

    const int n = Center + mesh_size[0];
    return field[n];
  }

  Vec3 GetZPlusFieldXZPeriodic(const int Z, const int Center) const noexcept
  {
    const int n = (Z + 1 == mesh_size[2])
                    ? Center - GetMeshSizeXZ() + (mesh_size[0] * 2)
                    : Center + mesh_size[0];
    return field[n];
  }

  Vec3 GetYMinusFieldYZ(const int Y, const int Center) const noexcept
  {
    if (Y == 0) return Vec3::Zero();

    const int n = Center - 1;
    return field[n];
  }

  Vec3 GetYMinusFieldYZPeriodic(const int Y, const int Center) const noexcept
  {
    const int n = (Y == 0) ? Center + mesh_size[1] - 2 : Center - 1;
    return field[n];
  }

  Vec3 GetYPlusFieldYZ(const int Y, const int Center) const noexcept
  {
    if (Y + 1 == mesh_size[1]) return Vec3::Zero();

    const int n = Center + 1;
    return field[n];
  }

  Vec3 GetYPlusFieldYZPeriodic(const int Y, const int Center) const noexcept
  {
    const int n =
      (Y + 1 == mesh_size[1]) ? Center - mesh_size[1] + 2 : Center + 1;
    return field[n];
  }

  Vec3 GetZMinusFieldYZ(const int Z, const int Center) const noexcept
  {
    if (Z == 0) return Vec3::Zero();

    const int n = Center - mesh_size[1];
    return field[n];
  }

  Vec3 GetZMinusFieldYZPeriodic(const int Z, const int Center) const noexcept
  {
    const int n = (Z == 0) ? Center + GetMeshSizeYZ() - (mesh_size[1] * 2)
                           : Center - mesh_size[1];
    return field[n];
  }

  Vec3 GetZPlusFieldYZ(const int Z, const int Center) const noexcept
  {
    if (Z + 1 == mesh_size[2]) return Vec3::Zero();

    const int n = Center + mesh_size[1];
    return field[n];
  }

  Vec3 GetZPlusFieldYZPeriodic(const int Z, const int Center) const noexcept
  {
    const int n = (Z + 1 == mesh_size[2])
                    ? Center - GetMeshSizeYZ() + (mesh_size[1] * 2)
                    : Center + mesh_size[1];
    return field[n];
  }

  Vec3 GetZMinusField3D(const int Z, const int Center) const noexcept
  {
    if (Z == 0) return Vec3::Zero();

    const int n = Center - GetMeshSizeXY();
    return field[n];
  }

  Vec3 GetZMinusField3DPeriodic(const int Z, const int Center) const noexcept
  {
    const int n = (Z == 0) ? Center + GetMeshSizeXYZ() - (GetMeshSizeXY() * 2)
                           : Center - GetMeshSizeXY();
    return field[n];
  }

  Vec3 GetZPlusField3D(const int Z, const int Center) const noexcept
  {
    if (Z + 1 == mesh_size[2]) return Vec3::Zero();

    const int n = Center + GetMeshSizeXY();
    return field[n];
  }

  Vec3 GetZPlusField3DPeriodic(const int Z, const int Center) const noexcept
  {
    const int n = (Z + 1 == mesh_size[2])
                    ? Center - GetMeshSizeXYZ() + (GetMeshSizeXY() * 2)
                    : Center + GetMeshSizeXY();
    return field[n];
  }

  bool isPtclInField2DXY(const std::array<int, 2>& XY) const noexcept
  {
    return ((0 < XY[0]) && (XY[0] < mesh_size[0]) && (0 < XY[1]) &&
            (XY[1] < mesh_size[1]))
             ? true
             : false;
  }

  int GetTgtMesh2DXY(const std::array<int, 2>& XY) const noexcept
  {
    return (mesh_size[0] * XY[1]) + XY[0];
  }

  bool isPtclInField2DXZ(const std::array<int, 2>& XZ) const noexcept
  {
    return ((0 < XZ[0]) && (XZ[0] < mesh_size[0]) && (0 < XZ[1]) &&
            (XZ[1] < mesh_size[2]))
             ? true
             : false;
  }

  int GetTgtMesh2DXZ(const std::array<int, 2>& XZ) const noexcept
  {
    return (mesh_size[0] * XZ[1]) + XZ[0];
  }

  bool isPtclInField2DYZ(const std::array<int, 2>& YZ) const noexcept
  {
    return ((0 < YZ[0]) && (YZ[0] < mesh_size[1]) && (0 < YZ[1]) &&
            (YZ[1] < mesh_size[2]))
             ? true
             : false;
  }

  int GetTgtMesh2DYZ(const std::array<int, 2>& YZ) const noexcept
  {
    return (mesh_size[1] * YZ[1]) + YZ[0];
  }

  bool isPtclInField3D(const std::array<int, 3>& XYZ) const noexcept
  {
    return ((0 < XYZ[0]) && (XYZ[0] < mesh_size[0]) && (0 < XYZ[1]) &&
            (XYZ[1] < mesh_size[1]) && (0 < XYZ[2]) && (XYZ[2] < mesh_size[2]))
             ? true
             : false;
  }

  int GetTgtMesh3D(const std::array<int, 3>& XYZ) const noexcept
  {
    return (GetMeshSizeXY() * XYZ[2]) + (mesh_size[0] * XYZ[1]) + XYZ[0];
  }
};

}  // namespace dem
