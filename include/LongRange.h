#pragma once
#include <cmath>
#include <vector>

#include "Field.h"
#include "LongRangeParticleList.h"
#include "Object.h"
#include "Particle.h"
#include "Vector.h"

namespace dem
{
class LongRange
{
private:
  Particle* ptcl;
  Object* obj;
  Field* field;
  LongRangeParticleList* lr_ptcl_list;

  double sum_energy;
  double pre_sum_energy;
  constexpr static double PI = 3.14159;
  constexpr static double EPS_E = 1.0e-18;
  constexpr static double EPS_B = 1.0e-18;

public:
  LongRange() noexcept
    : ptcl(nullptr),
      obj(nullptr),
      field(nullptr),
      lr_ptcl_list(nullptr),
      sum_energy(0.0),
      pre_sum_energy(0.0)
  {
  }

  void AssociateObj(Object* Obj) noexcept { obj = Obj; }

  void AssociateField(Field* Field) noexcept { field = Field; }

  void AssociatePtcl(Particle* Ptcl) noexcept { ptcl = Ptcl; }

  void AssociateLRPtcl(LongRangeParticleList* Long_Range_Ptcl) noexcept
  {
    lr_ptcl_list = Long_Range_Ptcl;
  }

  void CalcEleDipoleInteract(const double Dipole_Coef,
                             const double Field_Coef) noexcept;

  void CalcMagDipoleInteract(const double Dipole_Coef,
                             const double Field_Coef) noexcept;

  void SetFieldFromOtherDipole(const double Field_Coef) noexcept;
};

}  // namespace dem
