#include "LongRangeParticleList.h"

namespace dem
{
void LongRangeParticleList::ResizeArraySize() noexcept
{
  lr_ptcl_info.resize(ptcl->GetPtclStaySize());
}

void LongRangeParticleList::ModifyPtcl() noexcept
{
  if (ptcl->GetPtclStaySize() < GetLRPtclListSize())
  {
    lr_ptcl_info.clear();
    ResizeArraySize();
  }
}

void LongRangeParticleList::SetPtclInfo() noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

  std::vector<NeighbPtclInfo> ptcl_list;
  const Vec3 cell_len = lr_cell->GetCellLength() / ptcl->GetMinDiameter();
  const int max_num = static_cast<int>(cell_len.x * cell_len.y * cell_len.z);
  ptcl_list.reserve(max_num);

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    lr_ptcl_info[i].lr_ptcl.clear();
    const Vec3 p_post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();
    lr_cell->GetNeighbPtclList(i, ptcl_list);

    for (const auto& v : ptcl_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - p_post_pos_A + offset;
      const double post_dist = r_post_pos.Length();
      const Vec3 u_n_vec = r_post_pos / post_dist;
      const double R_sum = R_A + p[id_B].GetRadius();

      const double post_dist_corr = (post_dist <= R_sum) ? R_sum : post_dist;
      const Vec3 rel_post_pos_corr = post_dist_corr * u_n_vec;
      const double post_dist_inv = 1.0 / post_dist_corr;
      const double post_dist_inv_sq = post_dist_inv * post_dist_inv;
      const double d_inv_3 = post_dist_inv * post_dist_inv_sq;
      const double d_inv_5 = d_inv_3 * post_dist_inv_sq;
      const double d_inv_7 = d_inv_5 * post_dist_inv_sq;

#pragma omp critical
      {
        lr_ptcl_info[i].lr_ptcl.emplace_back(id_B, d_inv_3, d_inv_5, d_inv_7,
                                             -rel_post_pos_corr);

        lr_ptcl_info[id_B].lr_ptcl.emplace_back(i, d_inv_3, d_inv_5, d_inv_7,
                                                rel_post_pos_corr);
      }
    }
  }

#pragma omp single
  is_lr_ptcl_list_updated = true;
}

void LongRangeParticleList::SetPtclInfoCalcCoulombPtcl(
  const double Coef) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

  std::vector<NeighbPtclInfo> ptcl_list;
  const Vec3 cell_len = lr_cell->GetCellLength() / ptcl->GetMinDiameter();
  const int max_num = static_cast<int>(cell_len.x * cell_len.y * cell_len.z);
  ptcl_list.reserve(max_num);

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    lr_ptcl_info[i].lr_ptcl.clear();
    const Vec3 p_post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();
    lr_cell->GetNeighbPtclList(i, ptcl_list);

    for (const auto& v : ptcl_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - p_post_pos_A + offset;
      const double post_dist = r_post_pos.Length();
      const Vec3 u_n_vec = r_post_pos / post_dist;
      const double R_sum = R_A + p[id_B].GetRadius();

      const double post_dist_corr = (post_dist <= R_sum) ? R_sum : post_dist;
      const Vec3 rel_post_pos_corr = post_dist_corr * u_n_vec;
      const double post_dist_inv = 1.0 / post_dist_corr;
      const double post_dist_inv_sq = post_dist_inv * post_dist_inv;
      const double d_inv_3 = post_dist_inv * post_dist_inv_sq;
      const double d_inv_5 = d_inv_3 * post_dist_inv_sq;
      const double d_inv_7 = d_inv_5 * post_dist_inv_sq;

#pragma omp critical
      {
        lr_ptcl_info[i].lr_ptcl.emplace_back(id_B, d_inv_3, d_inv_5, d_inv_7,
                                             -rel_post_pos_corr);

        lr_ptcl_info[id_B].lr_ptcl.emplace_back(i, d_inv_3, d_inv_5, d_inv_7,
                                                rel_post_pos_corr);
      }

      const double chg = p[i].GetChg() * p[id_B].GetChg();
      const Vec3 C_Force = (Coef * post_dist_inv_sq * chg) * u_n_vec;
      p[i].AddCoulombPtclForce(-C_Force);
      p[id_B].AddCoulombPtclForce(C_Force);
    }
  }

#pragma omp single
  is_lr_ptcl_list_updated = true;
}

void LongRangeParticleList::CalcCoulombPtcl(const double Coef) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

  std::vector<NeighbPtclInfo> ptcl_list;
  const Vec3 cell_len = lr_cell->GetCellLength() / ptcl->GetMinDiameter();
  const int max_num = static_cast<int>(cell_len.x * cell_len.y * cell_len.z);
  ptcl_list.reserve(max_num);

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const Vec3 p_post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();
    lr_cell->GetNeighbPtclList(i, ptcl_list);

    for (const auto& v : ptcl_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - p_post_pos_A + offset;
      const double post_dist = r_post_pos.Length();
      const Vec3 u_n_vec = r_post_pos / post_dist;
      const double R_sum = R_A + p[id_B].GetRadius();
      const double post_dist_inv =
        (post_dist <= R_sum) ? 1.0 / R_sum : 1.0 / post_dist;

      const double chg = p[i].GetChg() * p[id_B].GetChg();
      const Vec3 C_Force =
        (Coef * post_dist_inv * post_dist_inv * chg) * u_n_vec;
      p[i].AddCoulombPtclForce(-C_Force);
      p[id_B].AddCoulombPtclForce(C_Force);
    }
  }
}

}  // namespace dem
