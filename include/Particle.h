#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "DataIO.h"
#include "Object.h"
#include "ParticleInfo.h"
#include "Scope.h"
#include "Vector.h"

namespace dem
{
class Particle
{
private:
  Object* obj;
  DataIO data_io;
  Scope scope;

  bool is_ptcl_stay_size_changed;
  double min_diameter, max_diameter;
  double trans_energy, rot_energy;
  std::vector<ParticleInfo> ptcl_stay;
  std::vector<ParticleInfo> ptcl_out;
  constexpr static double PI = 3.14159;

public:
  Particle() noexcept
    : obj(nullptr),
      is_ptcl_stay_size_changed(false),
      min_diameter(0.0),
      max_diameter(0.0),
      trans_energy(0.0),
      rot_energy(0.0)
  {
  }

  std::vector<ParticleInfo>& GetPtclStay() noexcept { return ptcl_stay; }

  std::vector<ParticleInfo>& GetPtclOut() noexcept { return ptcl_out; }

  void AssociateObj(Object* Obj) noexcept { obj = Obj; }

  bool isPtclStaySizeChanged() const noexcept
  {
    return is_ptcl_stay_size_changed;
  }

  int GetPtclStaySize() const noexcept
  {
    return static_cast<int>(ptcl_stay.size());
  }

  int GetPtclOutSize() const noexcept
  {
    return static_cast<int>(ptcl_out.size());
  }

  int GetPtclTotalSize() const noexcept
  {
    return GetPtclStaySize() + GetPtclOutSize();
  }

  double GetMinDiameter() const noexcept { return min_diameter; }

  double GetMaxDiameter() const noexcept { return max_diameter; }

  void PtclStaySizeisNotChanged() noexcept
  {
    is_ptcl_stay_size_changed = false;
  }

  void SetMinDiameter(const double v) noexcept { min_diameter = v; }

  void SetMaxDiameter(const double v) noexcept { max_diameter = v; }

  void SortPtclStay() noexcept;

  void SortPtclOut() noexcept;

  void LoadPtcl(const std::string& File) noexcept;

  void CheckBoundaryX(const double DomainX) noexcept;

  void CheckBoundaryY(const double DomainY) noexcept;

  void CheckBoundaryZ(const double DomainZ) noexcept;

  void CheckBoundaryXX(const double DomainX) noexcept;

  void CheckBoundaryYY(const double DomainY) noexcept;

  void CheckBoundaryZZ(const double DomainZ) noexcept;

  void ClearAdhForce() noexcept;

  void ClearImageForce() noexcept;

  void ClearCoulombPtclForce() noexcept;

  void ClearCoulombExtForce() noexcept;

  void ClearDieleForce() noexcept;

  void ClearDieleTorque() noexcept;

  void ClearMagForce() noexcept;

  void ClearMagTorque() noexcept;

  void StorePreExtForce() noexcept;

  void AddGravityForce(const Vec3& Gravity) noexcept;

  void AddCoulombExtForce() noexcept;

  void AddCoulombPtclForce() noexcept;

  void AddDieleForce() noexcept;

  void AddImageForce() noexcept;

  void AddMagForce() noexcept;

  void AddAdhForce() noexcept;

  void AddAirDrag(AirDragInfo& Info) noexcept;

  void AddAirDragCunningham(AirDragInfo& Info) noexcept;

  void AddAirDragReynolds(AirDragInfo& Info) noexcept;

  void StorePreExtTorque() noexcept;

  void AddDieleTorque() noexcept;

  void AddMagTorque() noexcept;

  void TimeIntegHVeloPosVerlet(const double Time_Step) noexcept;

  void TimeIntegVeloVerlet(const double Time_Step) noexcept;

  void TimeIntegHAngVeloAngVerlet(const double Time_Step) noexcept;

  void TimeIntegAngVeloVerlet(const double Time_Step) noexcept;

  std::array<double, 2> GetKineticEnergy() noexcept;
};

}  // namespace dem
