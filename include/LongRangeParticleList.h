#pragma once
#include <cmath>
#include <vector>

#include "Cell.h"
#include "Particle.h"
#include "Vector.h"

namespace dem
{
struct LongRangePtcl
{
  int id;
  double dist_inv_3, dist_inv_5, dist_inv_7;
  Vec3 rel_pos;

  LongRangePtcl() noexcept = default;

  LongRangePtcl(const int ID, const double D3, const double D5, const double D7,
                const Vec3& RP) noexcept
    : id(ID), dist_inv_3(D3), dist_inv_5(D5), dist_inv_7(D7), rel_pos(RP)
  {
  }
};

struct LongRangePtclInfo
{
  std::vector<LongRangePtcl> lr_ptcl;

  LongRangePtclInfo() noexcept = default;
};

class LongRangeParticleList
{
private:
  Particle* ptcl;
  Cell* lr_cell;
  bool is_lr_ptcl_list_updated;
  std::vector<LongRangePtclInfo> lr_ptcl_info;

public:
  LongRangeParticleList() noexcept
    : ptcl(nullptr), lr_cell(nullptr), is_lr_ptcl_list_updated(false)
  {
  }

  std::vector<LongRangePtclInfo>& GetLRPtclInfo() noexcept
  {
    return lr_ptcl_info;
  }

  void AssociatePtcl(Particle* Ptcl) noexcept { ptcl = Ptcl; }

  void AssociateLRCell(Cell* LR_Cell) noexcept { lr_cell = LR_Cell; }

  bool isLRPtclListUpdated() const noexcept { return is_lr_ptcl_list_updated; }

  int GetLRPtclListSize() const noexcept
  {
    return static_cast<int>(lr_ptcl_info.size());
  }

  void LRPtclListisNotUpdated() noexcept { is_lr_ptcl_list_updated = false; }

  void ResizeArraySize() noexcept;

  void ModifyPtcl() noexcept;

  void SetPtclInfo() noexcept;

  void SetPtclInfoCalcCoulombPtcl(const double Coef) noexcept;

  void CalcCoulombPtcl(const double Coef) noexcept;
};

}  // namespace dem
