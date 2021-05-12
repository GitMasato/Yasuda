#include "Collision.h"

namespace dem
{
void Collision::ClearCollOrder() noexcept { coll_info.clear(); }

void Collision::DetectCollPtcl(const std::array<bool, 3>& is_Periodic,
                               const Vec3& Domain, const bool is_Openmp,
                               const int Openmp_Size) noexcept
{
  const int ptcl_stay_size = ptcl->GetPtclStaySize();
  const int process_size = is_Openmp ? Openmp_Size : 1;
  const int reserve_size = ptcl_stay_size / process_size * 6;  // 6 directions

  if (is_Periodic[0] && is_Periodic[1] && is_Periodic[2])
  {
    DetectPtclCollPeriodicXYZ(Domain, reserve_size);
  }
  else if (is_Periodic[0] && is_Periodic[1])
  {
    DetectPtclCollPeriodicXY(Domain, reserve_size);
  }
  else if (is_Periodic[0] && is_Periodic[2])
  {
    DetectPtclCollPeriodicXZ(Domain, reserve_size);
  }
  else if (is_Periodic[1] && is_Periodic[2])
  {
    DetectPtclCollPeriodicYZ(Domain, reserve_size);
  }
  else if (is_Periodic[0])
  {
    DetectPtclCollPeriodicX(Domain, reserve_size);
  }
  else if (is_Periodic[1])
  {
    DetectPtclCollPeriodicY(Domain, reserve_size);
  }
  else if (is_Periodic[2])
  {
    DetectPtclCollPeriodicZ(Domain, reserve_size);
  }
  else
  {
    DetectPtclColl(reserve_size);
  }
}

void Collision::DetectPtclColl(const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      const Vec3 r_pos = p[id_B].GetPos() - pos_A;
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    // coll_info.reserve(coll_info.size() + coll.size());
    // coll_info.insert(coll_info.end(), coll.begin(), coll.end());

    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectPtclCollPeriodicX(const Vec3& Domain,
                                        const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());
  const Vec3 h_domain_sq = (Domain * 0.5) * (Domain * 0.5);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      Vec3 r_pos = p[id_B].GetPos() - pos_A;
      SetPeriodicEffect(r_pos.x, Domain.x, h_domain_sq.x);
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectPtclCollPeriodicY(const Vec3& Domain,
                                        const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());
  const Vec3 h_domain_sq = (Domain * 0.5) * (Domain * 0.5);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      Vec3 r_pos = p[id_B].GetPos() - pos_A;
      SetPeriodicEffect(r_pos.y, Domain.y, h_domain_sq.y);
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectPtclCollPeriodicZ(const Vec3& Domain,
                                        const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());
  const Vec3 h_domain_sq = (Domain * 0.5) * (Domain * 0.5);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      Vec3 r_pos = p[id_B].GetPos() - pos_A;
      SetPeriodicEffect(r_pos.z, Domain.z, h_domain_sq.z);
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectPtclCollPeriodicXY(const Vec3& Domain,
                                         const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());
  const Vec3 h_domain_sq = (Domain * 0.5) * (Domain * 0.5);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      Vec3 r_pos = p[id_B].GetPos() - pos_A;
      SetPeriodicEffect(r_pos.x, Domain.x, h_domain_sq.x);
      SetPeriodicEffect(r_pos.y, Domain.y, h_domain_sq.y);
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectPtclCollPeriodicXZ(const Vec3& Domain,
                                         const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());
  const Vec3 h_domain_sq = (Domain * 0.5) * (Domain * 0.5);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      Vec3 r_pos = p[id_B].GetPos() - pos_A;
      SetPeriodicEffect(r_pos.x, Domain.x, h_domain_sq.x);
      SetPeriodicEffect(r_pos.z, Domain.z, h_domain_sq.z);
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectPtclCollPeriodicYZ(const Vec3& Domain,
                                         const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());
  const Vec3 h_domain_sq = (Domain * 0.5) * (Domain * 0.5);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      Vec3 r_pos = p[id_B].GetPos() - pos_A;
      SetPeriodicEffect(r_pos.y, Domain.y, h_domain_sq.y);
      SetPeriodicEffect(r_pos.z, Domain.z, h_domain_sq.z);
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectPtclCollPeriodicXYZ(const Vec3& Domain,
                                          const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<NeighbPtclInfo> potential_coll_list;
  potential_coll_list.reserve(GetMaxPotentialCollNum());
  const Vec3 h_domain_sq = (Domain * 0.5) * (Domain * 0.5);

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    cell->GetNeighbPtclList(i, potential_coll_list);
    const Vec3 h_velo_A = p[i].GetHStepVelo();
    const Vec3 pos_A = p[i].GetPos();
    const Vec3 post_pos_A = p[i].GetPostPos();
    const double R_A = p[i].GetRadius();

    for (const auto& v : potential_coll_list)
    {
      const auto [id_B, offset] = v;
      const Vec3 r_post_pos = p[id_B].GetPostPos() - post_pos_A + offset;
      const double post_dist_sq = r_post_pos.LengthSq();
      const double R_sum = R_A + p[id_B].GetRadius();
      const double R_sum_sq = R_sum * R_sum;

      if (post_dist_sq > R_sum_sq) continue;  // no collision

      p[i].AddCollNumPtcl();
      p[id_B].AddCollNumPtcl();
      Vec3 r_pos = p[id_B].GetPos() - pos_A;
      SetPeriodicEffect(r_pos.x, Domain.x, h_domain_sq.x);
      SetPeriodicEffect(r_pos.y, Domain.y, h_domain_sq.y);
      SetPeriodicEffect(r_pos.z, Domain.z, h_domain_sq.z);
      const double dist_sq_diff = r_pos.LengthSq() - R_sum_sq;
      const auto [overlap, u_n_vec] =
        GetOverlapVec(post_dist_sq, R_sum, r_post_pos);

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(0, i, id_B, 0.0, overlap, u_n_vec, Vec3::Zero());
      }
      else
      {
        const Vec3 r_h_velo = p[id_B].GetHStepVelo() - h_velo_A;
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(0, i, id_B, time, overlap, u_n_vec, Vec3::Zero());
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectCollObj(const bool is_Openmp,
                              const int Openmp_Size) noexcept
{
  const int ptcl_stay_size = ptcl->GetPtclStaySize();
  const int process_size = is_Openmp ? Openmp_Size : 1;
  const int reserve_size = ptcl_stay_size / process_size;

  DetectFloorColl(reserve_size);
  DetectBoxColl(reserve_size);
  DetectHollowBoxColl(reserve_size);
  DetectCylinderColl(reserve_size);
  DetectHollowCylinderColl(reserve_size);
  DetectSphereColl(reserve_size);
  DetectHollowSphereColl(reserve_size);
}

void Collision::DetectCollObjImage(const bool is_Openmp, const int Openmp_Size,
                                   const double Image_Coef) noexcept
{
  const int ptcl_stay_size = ptcl->GetPtclStaySize();
  const int process_size = is_Openmp ? Openmp_Size : 1;
  const int reserve_size = ptcl_stay_size / process_size;

  DetectFloorCollImage(Image_Coef, reserve_size);
  DetectBoxCollImage(Image_Coef, reserve_size);
  DetectHollowBoxCollImage(Image_Coef, reserve_size);
  DetectCylinderCollImage(Image_Coef, reserve_size);
  DetectHollowCylinderCollImage(Image_Coef, reserve_size);
  DetectSphereCollImage(Image_Coef, reserve_size);
  DetectHollowSphereCollImage(Image_Coef, reserve_size);
}

void Collision::DetectFloorColl(const int Reserve) noexcept
{
  if (!(obj->GetFloorBaseSize())) return;

  std::vector<ObjectInfo>& o = obj->GetFloorBase();
  if (!(o[0].isCollOn())) return;

  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);
  constexpr Vec3 u_n_vec = -Vec3::UnitZ();

  const Vec3 o_h_velo = o[0].GetHStepVelo();
  const double o_pos_z = o[0].GetPosZ();
  const double o_post_pos_z = o[0].GetPostPosZ();

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for nowait schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const double p_R = p[i].GetRadius();
    const double post_surf_dist = p[i].GetPostPosZ() - p_R - o_post_pos_z;

    if (post_surf_dist > 0.0) continue;  // no collision

    p[i].AddCollNumObj();
    const double surf_dist = p[i].GetPosZ() - p_R - o_pos_z;

    if (surf_dist <= 0.0)
    {
      coll.emplace_back(1, i, 0, 0.0, -post_surf_dist, u_n_vec, o_h_velo);
    }
    else
    {
      const double r_h_velo = o_h_velo.z - p[i].GetHStepVeloZ();
      const double time = surf_dist / r_h_velo;
      coll.emplace_back(1, i, 0, time, -post_surf_dist, u_n_vec, o_h_velo);
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectFloorCollImage(const double Image_Coef,
                                     const int Reserve) noexcept
{
  if (!(obj->GetFloorBaseSize())) return;

  std::vector<ObjectInfo>& o = obj->GetFloorBase();
  if (!(o[0].isCollOn())) return;

  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);
  constexpr Vec3 u_n_vec = -Vec3::UnitZ();

  const Vec3 o_h_velo = o[0].GetHStepVelo();
  const double o_pos_z = o[0].GetPosZ();
  const double o_post_pos_z = o[0].GetPostPosZ();

  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for nowait schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const double p_R = p[i].GetRadius();
    const double r_post_pos = p[i].GetPostPosZ() - o_post_pos_z;
    const double post_surf_dist = r_post_pos - p_R;

    if (post_surf_dist > 0.0)  // no collision
    {
      p[i].AddImageForceZ(
        -GetImageForceLen(r_post_pos, Image_Coef, p[i].GetChgSq()));
      continue;
    }

    p[i].AddImageForceZ(-GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq()));
    p[i].AddCollNumObj();
    const double surf_dist = p[i].GetPosZ() - p_R - o_pos_z;

    if (surf_dist <= 0.0)
    {
      coll.emplace_back(1, i, 0, 0.0, -post_surf_dist, u_n_vec, o_h_velo);
    }
    else
    {
      const double r_h_velo = o_h_velo.z - p[i].GetHStepVeloZ();
      const double time = surf_dist / r_h_velo;
      coll.emplace_back(1, i, 0, time, -post_surf_dist, u_n_vec, o_h_velo);
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectBoxColl(const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<BoxInfo>& o = obj->GetBox();
  const int box_size = obj->GetBoxSize();

  for (int k = 0; k < box_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_pos_min = o[k].GetPos();
    const Vec3 o_pos_max = o_pos_min + o[k].GetDim();
    const Vec3 o_post_pos_min = o[k].GetPostPos();
    const Vec3 o_post_pos_max = o_post_pos_min + o[k].GetDim();

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 p_post_pos = p[i].GetPostPos();
      Vec3 o_post_coll_pt = p_post_pos;
      const bool is_in_box =
        ClampPos(o_post_coll_pt, o_post_pos_min, o_post_pos_max);

      // detect collision
      if (is_in_box)
      {
        p[i].AddCollNumObj();
        const auto [dir, post_dist] =
          GetCloseDirDist(p_post_pos, o_post_pos_min, o_post_pos_max);
        const double p_R = p[i].GetRadius();
        const double overlap = post_dist + p_R;
        const Vec3 u_n_vec = GetUNVecCloseDir(dir);
        const Vec3 p_pos = p[i].GetPos();
        const double dist = -GetDiffCloseDir(dir, p_pos, o_pos_min, o_pos_max);
        const double surf_dist = dist - p_R;

        if (surf_dist <= 0.0)
        {
          coll.emplace_back(2, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 p_h_velo = p[i].GetHStepVelo();
          const double r_h_velo =
            GetDiffCloseDir(dir, p_h_velo, o_h_velo, p_h_velo);
          const double time = surf_dist / r_h_velo;
          coll.emplace_back(2, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
      else
      {
        const Vec3 r_post_pos = o_post_coll_pt - p_post_pos;
        const double post_dist_sq = r_post_pos.LengthSq();
        const double p_R = p[i].GetRadius();
        const double p_R_sq = p_R * p_R;

        if (post_dist_sq > p_R_sq) continue;  // no collision

        p[i].AddCollNumObj();
        const auto [overlap, u_n_vec] =
          GetOverlapVec(post_dist_sq, p_R, r_post_pos);
        const Vec3 o_coll_pt = o_pos_min + o_post_coll_pt - o_post_pos_min;
        const Vec3 r_pos = o_coll_pt - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(2, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(2, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectBoxCollImage(const double Image_Coef,
                                   const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<BoxInfo>& o = obj->GetBox();
  const int box_size = obj->GetBoxSize();

  for (int k = 0; k < box_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_pos_min = o[k].GetPos();
    const Vec3 o_pos_max = o_pos_min + o[k].GetDim();
    const Vec3 o_post_pos_min = o[k].GetPostPos();
    const Vec3 o_post_pos_max = o_post_pos_min + o[k].GetDim();

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 p_post_pos = p[i].GetPostPos();
      Vec3 o_post_coll_pt = p_post_pos;
      const bool is_in_box =
        ClampPos(o_post_coll_pt, o_post_pos_min, o_post_pos_max);

      // detect collision
      if (is_in_box)
      {
        p[i].AddCollNumObj();
        const auto [dir, post_dist] =
          GetCloseDirDist(p_post_pos, o_post_pos_min, o_post_pos_max);
        const double p_R = p[i].GetRadius();
        const double overlap = post_dist + p_R;
        const Vec3 u_n_vec = GetUNVecCloseDir(dir);

        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const Vec3 p_pos = p[i].GetPos();
        const double dist = -GetDiffCloseDir(dir, p_pos, o_pos_min, o_pos_max);
        const double surf_dist = dist - p_R;

        if (surf_dist <= 0.0)
        {
          coll.emplace_back(2, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 p_h_velo = p[i].GetHStepVelo();
          const double r_h_velo =
            GetDiffCloseDir(dir, p_h_velo, o_h_velo, o_h_velo);
          const double time = surf_dist / r_h_velo;
          coll.emplace_back(2, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
      else
      {
        const Vec3 r_post_pos = o_post_coll_pt - p_post_pos;
        const double post_dist_sq = r_post_pos.LengthSq();
        const double p_R = p[i].GetRadius();
        const double p_R_sq = p_R * p_R;

        if (post_dist_sq > p_R_sq)  // no collision
        {
          const double I = Image_Coef * p[i].GetChgSq() / (4.0 * post_dist_sq);
          p[i].AddImageForce(r_post_pos / sqrt(post_dist_sq) * I);
          continue;
        }

        p[i].AddCollNumObj();
        const auto [overlap, u_n_vec] =
          GetOverlapVec(post_dist_sq, p_R, r_post_pos);

        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const Vec3 o_coll_pt = o_pos_min + o_post_coll_pt - o_post_pos_min;
        const Vec3 r_pos = o_coll_pt - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(2, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(2, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectHollowBoxColl(const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<HollowBoxInfo>& o = obj->GetHollowBox();
  const int hollow_box_size = obj->GetHollowBoxSize();

  for (int k = 0; k < hollow_box_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const double o_thk = o[k].GetThickness();
    const double o_thk_half = o_thk * 0.5;
    const Vec3 o_pos_min = o[k].GetPos();
    const Vec3 o_pos_max = o_pos_min + o[k].GetDim();
    const Vec3 o_pos_min_in = o_pos_min + o_thk;
    const Vec3 o_pos_max_in = o_pos_max - o_thk;
    const Vec3 o_post_pos_min = o[k].GetPostPos();
    const Vec3 o_post_pos_max = o_post_pos_min + o[k].GetDim();

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 p_post_pos = p[i].GetPostPos();
      Vec3 o_post_coll_pt = p_post_pos;
      const bool is_in_box =
        ClampPos(o_post_coll_pt, o_post_pos_min, o_post_pos_max);

      if (is_in_box)
      {
        const double p_R = p[i].GetRadius();
        const auto [dir, post_dist] =
          GetCloseDirDist(p_post_pos, o_post_pos_min, o_post_pos_max);

        // detect collision with outside wall
        if (post_dist <= o_thk_half)
        {
          p[i].AddCollNumObj();
          const double overlap = post_dist + p_R;
          const Vec3 u_n_vec = GetUNVecCloseDir(dir);
          const Vec3 p_pos = p[i].GetPos();
          const double dist =
            -GetDiffCloseDir(dir, p_pos, o_pos_min, o_pos_max);
          const double surf_dist = dist - p_R;

          if (surf_dist <= 0.0)
          {
            coll.emplace_back(3, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec3 p_h_velo = p[i].GetHStepVelo();
            const double r_h_velo =
              GetDiffCloseDir(dir, p_h_velo, o_h_velo, o_h_velo);
            const double time = surf_dist / r_h_velo;
            coll.emplace_back(3, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
        // detect collision with inside wall
        else if (post_dist <= o_thk + p_R)
        {
          p[i].AddCollNumObj();
          const double overlap = o_thk + p_R - post_dist;
          const Vec3 u_n_vec = -GetUNVecCloseDir(dir);
          const Vec3 p_pos = p[i].GetPos();
          const double dist =
            GetDiffCloseDir(dir, p_pos, o_pos_min_in, o_pos_max_in);
          const double surf_dist = dist - p_R;

          if (surf_dist <= 0.0)
          {
            coll.emplace_back(3, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec3 p_h_velo = p[i].GetHStepVelo();
            const double r_h_velo =
              -GetDiffCloseDir(dir, p_h_velo, o_h_velo, o_h_velo);
            const double time = surf_dist / r_h_velo;
            coll.emplace_back(3, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
      }
      else
      {
        const Vec3 r_post_pos = o_post_coll_pt - p_post_pos;
        const double post_dist_sq = r_post_pos.LengthSq();
        const double p_R = p[i].GetRadius();
        const double p_R_sq = p_R * p_R;

        if (post_dist_sq > p_R_sq) continue;  // no collision

        p[i].AddCollNumObj();
        const auto [overlap, u_n_vec] =
          GetOverlapVec(post_dist_sq, p_R, r_post_pos);
        const Vec3 o_coll_pt = o_pos_min + o_post_coll_pt - o_post_pos_min;
        const Vec3 r_pos = o_coll_pt - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(3, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(3, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectHollowBoxCollImage(const double Image_Coef,
                                         const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<HollowBoxInfo>& o = obj->GetHollowBox();
  const int hollow_box_size = obj->GetHollowBoxSize();

  for (int k = 0; k < hollow_box_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const double o_thk = o[k].GetThickness();
    const double o_thk_half = o_thk * 0.5;
    const Vec3 o_pos_min = o[k].GetPos();
    const Vec3 o_pos_max = o_pos_min + o[k].GetDim();
    const Vec3 o_pos_min_in = o_pos_min + o_thk;
    const Vec3 o_pos_max_in = o_pos_max - o_thk;
    const Vec3 o_post_pos_min = o[k].GetPostPos();
    const Vec3 o_post_pos_max = o_post_pos_min + o[k].GetDim();

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 p_post_pos = p[i].GetPostPos();
      Vec3 o_post_coll_pt = p_post_pos;
      const bool is_in_box =
        ClampPos(o_post_coll_pt, o_post_pos_min, o_post_pos_max);

      if (is_in_box)
      {
        const double p_R = p[i].GetRadius();
        const auto [dir, post_dist] =
          GetCloseDirDist(p_post_pos, o_post_pos_min, o_post_pos_max);

        // detect collision with outside wall
        if (post_dist <= o_thk_half)
        {
          p[i].AddCollNumObj();
          const double overlap = post_dist + p_R;
          const Vec3 u_n_vec = GetUNVecCloseDir(dir);

          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec * I);

          const Vec3 p_pos = p[i].GetPos();
          const double dist =
            -GetDiffCloseDir(dir, p_pos, o_pos_min, o_pos_max);
          const double surf_dist = dist - p_R;

          if (surf_dist <= 0.0)
          {
            coll.emplace_back(3, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec3 p_h_velo = p[i].GetHStepVelo();
            const double r_h_velo =
              GetDiffCloseDir(dir, p_h_velo, o_h_velo, o_h_velo);
            const double time = surf_dist / r_h_velo;
            coll.emplace_back(3, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
        // detect collision with inside wall
        else if (post_dist <= o_thk + p_R)
        {
          p[i].AddCollNumObj();
          const double overlap = o_thk + p_R - post_dist;
          const Vec3 u_n_vec = -GetUNVecCloseDir(dir);

          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec * I);

          const Vec3 p_pos = p[i].GetPos();
          const double dist =
            GetDiffCloseDir(dir, p_pos, o_pos_min_in, o_pos_max_in);
          const double surf_dist = dist - p_R;

          if (surf_dist <= 0.0)
          {
            coll.emplace_back(3, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec3 p_h_velo = p[i].GetHStepVelo();
            const double r_h_velo =
              -GetDiffCloseDir(dir, p_h_velo, o_h_velo, o_h_velo);
            const double time = surf_dist / r_h_velo;
            coll.emplace_back(3, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
        // inside box, no collision
        else
        {
          const double I_dist = post_dist - o_thk;
          const double I =
            GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(-GetUNVecCloseDir(dir) * I);
        }
      }
      else
      {
        const Vec3 r_post_pos = o_post_coll_pt - p_post_pos;
        const double post_dist_sq = r_post_pos.LengthSq();
        const double p_R = p[i].GetRadius();
        const double p_R_sq = p_R * p_R;

        if (post_dist_sq > p_R_sq)  // no collision
        {
          const double I = Image_Coef * p[i].GetChgSq() / (4.0 * post_dist_sq);
          p[i].AddImageForce(r_post_pos / sqrt(post_dist_sq) * I);
          continue;
        }

        p[i].AddCollNumObj();
        const auto [overlap, u_n_vec] =
          GetOverlapVec(post_dist_sq, p_R, r_post_pos);

        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const Vec3 o_coll_pt = o_pos_min + o_post_coll_pt - o_post_pos_min;
        const Vec3 r_pos = o_coll_pt - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(3, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(3, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectCylinderColl(const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<CylinderInfo>& o = obj->GetCylinder();
  const int cylinder_size = obj->GetCylinderSize();

  for (int k = 0; k < cylinder_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_h_velo_rot = o[k].Rotated(o_h_velo);
    const Vec3 o_pos_rot = o[k].Rotated(o[k].GetPos());
    const Vec3 u_n_vec_x_corr = o[k].RotatedReturn(Vec3::UnitX());
    const Vec3 o_post_pos_rot = o[k].Rotated(o[k].GetPostPos());

    const double o_len = o[k].GetLength();
    const double o_post_pos_rot_x_mid = o_post_pos_rot.x + (o_len * 0.5);
    const double o_outer_R = o[k].GetOuterRadius();
    const double o_outer_R_sq = o_outer_R * o_outer_R;

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const double p_R = p[i].GetRadius();
      const Vec3 p_post_pos_rot = o[k].Rotated(p[i].GetPostPos());
      Vec3 r_post_pos_rot = o_post_pos_rot - p_post_pos_rot;
      const double post_dist_yz_rot_sq = r_post_pos_rot.YZ().LengthSq();

      if (post_dist_yz_rot_sq <= o_outer_R_sq)
      {
        if (p_post_pos_rot.x <= o_post_pos_rot_x_mid)
        {
          // no collision
          if (r_post_pos_rot.x > p_R) continue;

          p[i].AddCollNumObj();
          const double overlap = -r_post_pos_rot.x + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = o_pos_rot.x - p_pos_rot - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(4, i, k, 0.0, overlap, u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = p_h_velo_rot - o_h_velo_rot.x;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(4, i, k, time, overlap, u_n_vec_x_corr, o_h_velo);
          }
        }
        else
        {
          // no collision
          if (-(r_post_pos_rot.x + o_len) > p_R) continue;

          p[i].AddCollNumObj();
          const double overlap = r_post_pos_rot.x + o_len + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = p_pos_rot - (o_pos_rot.x + o_len) - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(4, i, k, 0.0, overlap, -u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = o_h_velo_rot.x - p_h_velo_rot;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(4, i, k, time, overlap, -u_n_vec_x_corr,
                              o_h_velo);
          }
        }
      }
      else if ((o_post_pos_rot.x <= p_post_pos_rot.x) &&
               (p_post_pos_rot.x <= o_post_pos_rot.x + o_len))
      {
        const double R_outer_sum = o_outer_R + p_R;
        const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

        // no collision
        if (post_dist_yz_rot_sq > R_outer_sum_sq) continue;

        p[i].AddCollNumObj();
        r_post_pos_rot.x = 0.0;
        const double post_dist = sqrt(post_dist_yz_rot_sq);
        const double overlap = R_outer_sum - post_dist;
        const Vec3 u_n_vec_rot = r_post_pos_rot / post_dist;
        const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

        const Vec2 p_pos_rot = o[k].Rotated(p[i].GetPos()).YZ();
        const Vec2 r_pos_rot = o_pos_rot.YZ() - p_pos_rot;
        const double dist_sq_diff = r_pos_rot.LengthSq() - R_outer_sum_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(4, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec2 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo()).YZ();
          const Vec2 r_h_velo_rot = o_h_velo_rot.YZ() - p_h_velo_rot;
          const double time =
            GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
          coll.emplace_back(4, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
      else
      {
        r_post_pos_rot.x = 0.0;
        Vec3 pos_corr_rot =
          -(o_outer_R) * (r_post_pos_rot / sqrt(post_dist_yz_rot_sq));
        if (o_post_pos_rot_x_mid <= p_post_pos_rot.x) pos_corr_rot.x += o_len;

        const Vec3 r_post_pos_coll_pt_rot =
          o_post_pos_rot + pos_corr_rot - p_post_pos_rot;
        const double post_dist_sq = r_post_pos_coll_pt_rot.LengthSq();
        const double p_R_sq = p_R * p_R;

        // no collision
        if (post_dist_sq > p_R_sq) continue;

        p[i].AddCollNumObj();
        const double post_dist = sqrt(post_dist_sq);
        const double overlap = p_R - post_dist;
        const Vec3 u_n_vec_rot = r_post_pos_coll_pt_rot / post_dist;
        const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

        const Vec3 p_pos_rot = o[k].Rotated(p[i].GetPos());
        const Vec3 r_pos_rot = o_pos_rot + pos_corr_rot - p_pos_rot;
        const double dist_sq_diff = r_pos_rot.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(4, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo());
          const Vec3 r_h_velo_rot = o_h_velo_rot - p_h_velo_rot;
          const double time =
            GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
          coll.emplace_back(4, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectCylinderCollImage(const double Image_Coef,
                                        const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<CylinderInfo>& o = obj->GetCylinder();
  const int cylinder_size = obj->GetCylinderSize();

  for (int k = 0; k < cylinder_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_h_velo_rot = o[k].Rotated(o_h_velo);
    const Vec3 o_pos_rot = o[k].Rotated(o[k].GetPos());
    const Vec3 u_n_vec_x_corr = o[k].RotatedReturn(Vec3::UnitX());
    const Vec3 o_post_pos_rot = o[k].Rotated(o[k].GetPostPos());

    const double o_len = o[k].GetLength();
    const double o_post_pos_rot_x_mid = o_post_pos_rot.x + (o_len * 0.5);
    const double o_outer_R = o[k].GetOuterRadius();
    const double o_outer_R_sq = o_outer_R * o_outer_R;

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const double p_R = p[i].GetRadius();
      const Vec3 p_post_pos_rot = o[k].Rotated(p[i].GetPostPos());
      Vec3 r_post_pos_rot = o_post_pos_rot - p_post_pos_rot;
      const double post_dist_yz_rot_sq = r_post_pos_rot.YZ().LengthSq();

      if (post_dist_yz_rot_sq <= o_outer_R_sq)
      {
        if (p_post_pos_rot.x <= o_post_pos_rot_x_mid)
        {
          // no collision
          if (r_post_pos_rot.x > p_R)
          {
            const double I_dist = r_post_pos_rot.x;
            const double I =
              GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
            p[i].AddImageForce(u_n_vec_x_corr * I);
            continue;
          }

          p[i].AddCollNumObj();
          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec_x_corr * I);

          const double overlap = -r_post_pos_rot.x + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = o_pos_rot.x - p_pos_rot - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(4, i, k, 0.0, overlap, u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = p_h_velo_rot - o_h_velo_rot.x;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(4, i, k, time, overlap, u_n_vec_x_corr, o_h_velo);
          }
        }
        else
        {
          // no collision
          if (-(r_post_pos_rot.x + o_len) > p_R)
          {
            const double I_dist = -(r_post_pos_rot.x + o_len);
            const double I =
              GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
            p[i].AddImageForce(-u_n_vec_x_corr * I);
            continue;
          }

          p[i].AddCollNumObj();
          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(-u_n_vec_x_corr * I);

          const double overlap = r_post_pos_rot.x + o_len + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = p_pos_rot - (o_pos_rot.x + o_len) - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(4, i, k, 0.0, overlap, -u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = o_h_velo_rot.x - p_h_velo_rot;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(4, i, k, time, overlap, -u_n_vec_x_corr,
                              o_h_velo);
          }
        }
      }
      else if ((o_post_pos_rot.x <= p_post_pos_rot.x) &&
               (p_post_pos_rot.x <= o_post_pos_rot.x + o_len))
      {
        r_post_pos_rot.x = 0.0;
        const double post_dist = sqrt(post_dist_yz_rot_sq);
        const Vec3 u_n_vec_rot = r_post_pos_rot / post_dist;
        const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

        const double R_outer_sum = o_outer_R + p_R;
        const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

        // no collision
        if (post_dist_yz_rot_sq > R_outer_sum_sq)
        {
          const double I_dist = post_dist - o_outer_R;
          const double I =
            GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec * I);
          continue;
        }

        p[i].AddCollNumObj();
        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const double overlap = R_outer_sum - post_dist;
        const Vec2 p_pos_rot = o[k].Rotated(p[i].GetPos()).YZ();
        const Vec2 r_pos_rot = o_pos_rot.YZ() - p_pos_rot;
        const double dist_sq_diff = r_pos_rot.LengthSq() - R_outer_sum_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(4, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec2 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo()).YZ();
          const Vec2 r_h_velo_rot = o_h_velo_rot.YZ() - p_h_velo_rot;
          const double time =
            GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
          coll.emplace_back(4, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
      else
      {
        r_post_pos_rot.x = 0.0;
        Vec3 pos_corr_rot =
          -(o_outer_R) * (r_post_pos_rot / sqrt(post_dist_yz_rot_sq));
        if (o_post_pos_rot_x_mid <= p_post_pos_rot.x) pos_corr_rot.x += o_len;

        const Vec3 r_post_pos_coll_pt_rot =
          o_post_pos_rot + pos_corr_rot - p_post_pos_rot;
        const double post_dist_sq = r_post_pos_coll_pt_rot.LengthSq();
        const double p_R_sq = p_R * p_R;

        const double post_dist = sqrt(post_dist_sq);
        const Vec3 u_n_vec_rot = r_post_pos_coll_pt_rot / post_dist;
        const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

        // no collision
        if (post_dist_sq > p_R_sq)
        {
          const double I = Image_Coef * p[i].GetChgSq() / (4.0 * post_dist_sq);
          p[i].AddImageForce(u_n_vec * I);
          continue;
        }

        p[i].AddCollNumObj();
        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const double overlap = p_R - post_dist;
        const Vec3 p_pos_rot = o[k].Rotated(p[i].GetPos());
        const Vec3 r_pos_rot = o_pos_rot + pos_corr_rot - p_pos_rot;
        const double dist_sq_diff = r_pos_rot.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(4, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo());
          const Vec3 r_h_velo_rot = o_h_velo_rot - p_h_velo_rot;
          const double time =
            GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
          coll.emplace_back(4, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectHollowCylinderColl(const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<HollowCylinderInfo>& o = obj->GetHollowCylinder();
  const int hollow_cylinder_size = obj->GetHollowCylinderSize();

  for (int k = 0; k < hollow_cylinder_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_h_velo_rot = o[k].Rotated(o_h_velo);
    const Vec3 o_pos_rot = o[k].Rotated(o[k].GetPos());
    const Vec3 u_n_vec_x_corr = o[k].RotatedReturn(Vec3::UnitX());
    const Vec3 o_post_pos_rot = o[k].Rotated(o[k].GetPostPos());

    const double o_len = o[k].GetLength();
    const double o_post_pos_rot_x_mid = o_post_pos_rot.x + (o_len * 0.5);
    const double o_inner_R = o[k].GetInnerRadius();
    const double o_outer_R = o[k].GetOuterRadius();
    const double o_middle_R_sq =
      (o_inner_R + o_outer_R) * (o_inner_R + o_outer_R) * 0.25;
    const double o_inner_R_sq = o_inner_R * o_inner_R;
    const double o_outer_R_sq = o_outer_R * o_outer_R;

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const double p_R = p[i].GetRadius();
      const Vec3 p_post_pos_rot = o[k].Rotated(p[i].GetPostPos());
      Vec3 r_post_pos_rot = o_post_pos_rot - p_post_pos_rot;
      const double post_dist_yz_rot_sq = r_post_pos_rot.YZ().LengthSq();

      if ((o_inner_R_sq <= post_dist_yz_rot_sq) &&
          (post_dist_yz_rot_sq <= o_outer_R_sq))
      {
        if (p_post_pos_rot.x <= o_post_pos_rot_x_mid)
        {
          // no collision
          if (r_post_pos_rot.x > p_R) continue;

          p[i].AddCollNumObj();
          const double overlap = -r_post_pos_rot.x + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = o_pos_rot.x - p_pos_rot - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = p_h_velo_rot - o_h_velo_rot.x;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(5, i, k, time, overlap, u_n_vec_x_corr, o_h_velo);
          }
        }
        else
        {
          // no collision
          if (-(r_post_pos_rot.x + o_len) > p_R) continue;

          p[i].AddCollNumObj();
          const double overlap = r_post_pos_rot.x + o_len + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = p_pos_rot - (o_pos_rot.x + o_len) - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, -u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = o_h_velo_rot.x - p_h_velo_rot;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(5, i, k, time, overlap, -u_n_vec_x_corr,
                              o_h_velo);
          }
        }
      }
      else if ((o_post_pos_rot.x <= p_post_pos_rot.x) &&
               (p_post_pos_rot.x <= o_post_pos_rot.x + o_len))
      {
        if (post_dist_yz_rot_sq <= o_middle_R_sq)
        {
          const double R_inner_red = o_inner_R - p_R;
          const double R_inner_red_sq = R_inner_red * R_inner_red;

          // no collision
          if (post_dist_yz_rot_sq < R_inner_red_sq) continue;

          p[i].AddCollNumObj();
          r_post_pos_rot.x = 0.0;
          const double post_dist = sqrt(post_dist_yz_rot_sq);
          const double overlap = post_dist - R_inner_red;
          const Vec3 u_n_vec_rot = r_post_pos_rot / post_dist;
          const Vec3 u_n_vec = -o[k].RotatedReturn(u_n_vec_rot);

          const Vec2 p_pos_rot = o[k].Rotated(p[i].GetPos()).YZ();
          const Vec2 r_pos_rot = o_pos_rot.YZ() - p_pos_rot;
          const double dist_sq_diff = r_pos_rot.LengthSq() - R_inner_red_sq;

          if (dist_sq_diff >= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec2 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo()).YZ();
            const Vec2 r_h_velo_rot = o_h_velo_rot.YZ() - p_h_velo_rot;
            const double time =
              GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
            coll.emplace_back(5, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
        else
        {
          const double R_outer_sum = o_outer_R + p_R;
          const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

          // no collision
          if (R_outer_sum_sq < post_dist_yz_rot_sq) continue;

          p[i].AddCollNumObj();
          r_post_pos_rot.x = 0.0;
          const double post_dist = sqrt(post_dist_yz_rot_sq);
          const double overlap = R_outer_sum - post_dist;
          const Vec3 u_n_vec_rot = r_post_pos_rot / post_dist;
          const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

          const Vec2 p_pos_rot = o[k].Rotated(p[i].GetPos()).YZ();
          const Vec2 r_pos_rot = o_pos_rot.YZ() - p_pos_rot;
          const double dist_sq_diff = r_pos_rot.LengthSq() - R_outer_sum_sq;

          if (dist_sq_diff <= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec2 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo()).YZ();
            const Vec2 r_h_velo_rot = o_h_velo_rot.YZ() - p_h_velo_rot;
            const double time =
              GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
            coll.emplace_back(5, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
      }
      else
      {
        r_post_pos_rot.x = 0.0;
        Vec3 pos_corr_rot =
          post_dist_yz_rot_sq <= o_middle_R_sq
            ? -o_inner_R * r_post_pos_rot / sqrt(post_dist_yz_rot_sq)
            : -o_outer_R * r_post_pos_rot / sqrt(post_dist_yz_rot_sq);
        if (o_post_pos_rot_x_mid <= p_post_pos_rot.x) pos_corr_rot.x += o_len;

        const Vec3 r_post_pos_coll_pt_rot =
          o_post_pos_rot + pos_corr_rot - p_post_pos_rot;
        const double post_dist_sq = r_post_pos_coll_pt_rot.LengthSq();
        const double p_R_sq = p_R * p_R;

        // no collision
        if (post_dist_sq > p_R_sq) continue;

        p[i].AddCollNumObj();
        const double post_dist = sqrt(post_dist_sq);
        const double overlap = p_R - post_dist;
        const Vec3 u_n_vec_rot = r_post_pos_coll_pt_rot / post_dist;
        const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

        const Vec3 r_pos_rot =
          o_pos_rot + pos_corr_rot - o[k].Rotated(p[i].GetPos());
        const double dist_sq_diff = r_pos_rot.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo());
          const Vec3 r_h_velo_rot = o_h_velo_rot - p_h_velo_rot;
          const double time =
            GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
          coll.emplace_back(5, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectHollowCylinderCollImage(const double Image_Coef,
                                              const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<HollowCylinderInfo>& o = obj->GetHollowCylinder();
  const int hollow_cylinder_size = obj->GetHollowCylinderSize();

  for (int k = 0; k < hollow_cylinder_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_h_velo_rot = o[k].Rotated(o_h_velo);
    const Vec3 o_pos_rot = o[k].Rotated(o[k].GetPos());
    const Vec3 u_n_vec_x_corr = o[k].RotatedReturn(Vec3::UnitX());
    const Vec3 o_post_pos_rot = o[k].Rotated(o[k].GetPostPos());

    const double o_len = o[k].GetLength();
    const double o_post_pos_rot_x_mid = o_post_pos_rot.x + (o_len * 0.5);
    const double o_inner_R = o[k].GetInnerRadius();
    const double o_outer_R = o[k].GetOuterRadius();
    const double o_middle_R_sq =
      (o_inner_R + o_outer_R) * (o_inner_R + o_outer_R) * 0.25;
    const double o_inner_R_sq = o_inner_R * o_inner_R;
    const double o_outer_R_sq = o_outer_R * o_outer_R;

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const double p_R = p[i].GetRadius();
      const Vec3 p_post_pos_rot = o[k].Rotated(p[i].GetPostPos());
      Vec3 r_post_pos_rot = o_post_pos_rot - p_post_pos_rot;
      const double post_dist_yz_rot_sq = r_post_pos_rot.YZ().LengthSq();

      if ((o_inner_R_sq <= post_dist_yz_rot_sq) &&
          (post_dist_yz_rot_sq <= o_outer_R_sq))
      {
        if (p_post_pos_rot.x <= o_post_pos_rot_x_mid)
        {
          // no collision
          if (r_post_pos_rot.x > p_R)
          {
            const double I_dist = r_post_pos_rot.x;
            const double I =
              GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
            p[i].AddImageForce(u_n_vec_x_corr * I);
            continue;
          }

          p[i].AddCollNumObj();
          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec_x_corr * I);

          const double overlap = -r_post_pos_rot.x + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = o_pos_rot.x - p_pos_rot - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = p_h_velo_rot - o_h_velo_rot.x;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(5, i, k, time, overlap, u_n_vec_x_corr, o_h_velo);
          }
        }
        else
        {
          // no collision
          if (-(r_post_pos_rot.x + o_len) > p_R)
          {
            const double I_dist = -(r_post_pos_rot.x + o_len);
            const double I =
              GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
            p[i].AddImageForce(-u_n_vec_x_corr * I);
            continue;
          }

          p[i].AddCollNumObj();
          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(-u_n_vec_x_corr * I);

          const double overlap = r_post_pos_rot.x + o_len + p_R;
          const double p_pos_rot = o[k].RotatedX(p[i].GetPos());
          const double surf_dist_rot = p_pos_rot - (o_pos_rot.x + o_len) - p_R;

          if (surf_dist_rot <= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, -u_n_vec_x_corr, o_h_velo);
          }
          else
          {
            const double p_h_velo_rot = o[k].RotatedX(p[i].GetHStepVelo());
            const double r_h_velo_rot = o_h_velo_rot.x - p_h_velo_rot;
            const double time = surf_dist_rot / r_h_velo_rot;
            coll.emplace_back(5, i, k, time, overlap, -u_n_vec_x_corr,
                              o_h_velo);
          }
        }
      }
      else if ((o_post_pos_rot.x <= p_post_pos_rot.x) &&
               (p_post_pos_rot.x <= o_post_pos_rot.x + o_len))
      {
        if (post_dist_yz_rot_sq <= o_middle_R_sq)
        {
          r_post_pos_rot.x = 0.0;
          const double post_dist = sqrt(post_dist_yz_rot_sq);
          const Vec3 u_n_vec_rot = r_post_pos_rot / post_dist;
          const Vec3 u_n_vec = -o[k].RotatedReturn(u_n_vec_rot);

          const double R_inner_red = o_inner_R - p_R;
          const double R_inner_red_sq = R_inner_red * R_inner_red;

          // no collision
          if (post_dist_yz_rot_sq < R_inner_red_sq)
          {
            const double I_dist = o_inner_R - post_dist;
            const double I =
              GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
            p[i].AddImageForce(u_n_vec * I);
            continue;
          }

          p[i].AddCollNumObj();
          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec * I);

          const double overlap = post_dist - R_inner_red;
          const Vec2 p_pos_rot = o[k].Rotated(p[i].GetPos()).YZ();
          const Vec2 r_pos_rot = o_pos_rot.YZ() - p_pos_rot;
          const double dist_sq_diff = r_pos_rot.LengthSq() - R_inner_red_sq;

          if (dist_sq_diff >= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec2 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo()).YZ();
            const Vec2 r_h_velo_rot = o_h_velo_rot.YZ() - p_h_velo_rot;
            const double time =
              GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
            coll.emplace_back(5, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
        else
        {
          r_post_pos_rot.x = 0.0;
          const double post_dist = sqrt(post_dist_yz_rot_sq);
          const Vec3 u_n_vec_rot = r_post_pos_rot / post_dist;
          const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

          const double R_outer_sum = o_outer_R + p_R;
          const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

          // no collision
          if (R_outer_sum_sq < post_dist_yz_rot_sq)
          {
            const double I_dist = post_dist - o_outer_R;
            const double I =
              GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
            p[i].AddImageForce(u_n_vec * I);
            continue;
          }

          p[i].AddCollNumObj();
          const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec * I);

          const double overlap = R_outer_sum - post_dist;
          const Vec2 p_pos_rot = o[k].Rotated(p[i].GetPos()).YZ();
          const Vec2 r_pos_rot = o_pos_rot.YZ() - p_pos_rot;
          const double dist_sq_diff = r_pos_rot.LengthSq() - R_outer_sum_sq;

          if (dist_sq_diff <= 0.0)
          {
            coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec, o_h_velo);
          }
          else
          {
            const Vec2 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo()).YZ();
            const Vec2 r_h_velo_rot = o_h_velo_rot.YZ() - p_h_velo_rot;
            const double time =
              GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
            coll.emplace_back(5, i, k, time, overlap, u_n_vec, o_h_velo);
          }
        }
      }
      else
      {
        r_post_pos_rot.x = 0.0;
        Vec3 pos_corr_rot =
          post_dist_yz_rot_sq <= o_middle_R_sq
            ? -o_inner_R * r_post_pos_rot / sqrt(post_dist_yz_rot_sq)
            : -o_outer_R * r_post_pos_rot / sqrt(post_dist_yz_rot_sq);
        if (o_post_pos_rot_x_mid <= p_post_pos_rot.x) pos_corr_rot.x += o_len;

        const Vec3 r_post_pos_coll_pt_rot =
          o_post_pos_rot + pos_corr_rot - p_post_pos_rot;
        const double post_dist_sq = r_post_pos_coll_pt_rot.LengthSq();
        const double p_R_sq = p_R * p_R;

        const double post_dist = sqrt(post_dist_sq);
        const Vec3 u_n_vec_rot = r_post_pos_coll_pt_rot / post_dist;
        const Vec3 u_n_vec = o[k].RotatedReturn(u_n_vec_rot);

        // no collision
        if (post_dist_sq > p_R_sq)
        {
          const double I = Image_Coef * p[i].GetChgSq() / (4.0 * post_dist_sq);
          p[i].AddImageForce(u_n_vec * I);
          continue;
        }

        p[i].AddCollNumObj();
        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const double overlap = p_R - post_dist;
        const Vec3 r_pos_rot =
          o_pos_rot + pos_corr_rot - o[k].Rotated(p[i].GetPos());
        const double dist_sq_diff = r_pos_rot.LengthSq() - p_R_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(5, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 p_h_velo_rot = o[k].Rotated(p[i].GetHStepVelo());
          const Vec3 r_h_velo_rot = o_h_velo_rot - p_h_velo_rot;
          const double time =
            GetCollTime(r_pos_rot, r_h_velo_rot, dist_sq_diff);
          coll.emplace_back(5, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectSphereColl(const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<SphereInfo>& o = obj->GetSphere();
  const int sphere_size = obj->GetSphereSize();

  for (int k = 0; k < sphere_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_pos = o[k].GetPos();
    const Vec3 o_post_pos = o[k].GetPostPos();
    const double o_outer_R = o[k].GetOuterRadius();

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 r_post_pos = o_post_pos - p[i].GetPostPos();
      const double post_dist_sq = r_post_pos.LengthSq();
      const double p_R = p[i].GetRadius();
      const double R_outer_sum = o_outer_R + p_R;
      const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

      if (post_dist_sq > R_outer_sum_sq) continue;  // no collision

      p[i].AddCollNumObj();
      const double post_dist = sqrt(post_dist_sq);
      const Vec3 u_n_vec = r_post_pos / post_dist;
      const double overlap = R_outer_sum - post_dist;

      const Vec3 r_pos = o_pos - p[i].GetPos();
      const double dist_sq_diff = r_pos.LengthSq() - R_outer_sum_sq;

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(6, i, k, 0.0, overlap, u_n_vec, o_h_velo);
      }
      else
      {
        const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(6, i, k, time, overlap, u_n_vec, o_h_velo);
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectSphereCollImage(const double Image_Coef,
                                      const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<SphereInfo>& o = obj->GetSphere();
  const int sphere_size = obj->GetSphereSize();

  for (int k = 0; k < sphere_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_pos = o[k].GetPos();
    const Vec3 o_post_pos = o[k].GetPostPos();
    const double o_outer_R = o[k].GetOuterRadius();

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 r_post_pos = o_post_pos - p[i].GetPostPos();
      const double post_dist_sq = r_post_pos.LengthSq();
      const double p_R = p[i].GetRadius();
      const double R_outer_sum = o_outer_R + p_R;
      const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

      const double post_dist = sqrt(post_dist_sq);
      const Vec3 u_n_vec = r_post_pos / post_dist;

      if (post_dist_sq > R_outer_sum_sq)  // no collision
      {
        const double I_dist = post_dist - o_outer_R;
        const double I = GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);
        continue;
      }

      p[i].AddCollNumObj();
      const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
      p[i].AddImageForce(u_n_vec * I);

      const double overlap = R_outer_sum - post_dist;
      const Vec3 r_pos = o_pos - p[i].GetPos();
      const double dist_sq_diff = r_pos.LengthSq() - R_outer_sum_sq;

      if (dist_sq_diff <= 0.0)
      {
        coll.emplace_back(6, i, k, 0.0, overlap, u_n_vec, o_h_velo);
      }
      else
      {
        const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
        const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
        coll.emplace_back(6, i, k, time, overlap, u_n_vec, o_h_velo);
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectHollowSphereColl(const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<HollowSphereInfo>& o = obj->GetHollowSphere();
  const int hollow_sphere_size = obj->GetHollowSphereSize();

  for (int k = 0; k < hollow_sphere_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_pos = o[k].GetPos();
    const Vec3 o_post_pos = o[k].GetPostPos();
    const double o_inner_R = o[k].GetInnerRadius();
    const double o_outer_R = o[k].GetOuterRadius();
    const double o_middle_R_sq =
      (o_inner_R + o_outer_R) * (o_inner_R + o_outer_R) * 0.25;

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 r_post_pos = o_post_pos - p[i].GetPostPos();
      const double post_dist_sq = r_post_pos.LengthSq();
      const double p_R = p[i].GetRadius();

      if (post_dist_sq <= o_middle_R_sq)
      {
        const double R_inner_red = o_inner_R - p_R;
        const double R_inner_red_sq = R_inner_red * R_inner_red;

        if (R_inner_red_sq > post_dist_sq) continue;  // no collision

        p[i].AddCollNumObj();
        const double post_dist = sqrt(post_dist_sq);
        const Vec3 u_n_vec = -r_post_pos / post_dist;
        const double overlap = post_dist - R_inner_red;

        const Vec3 r_pos = o_pos - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - R_inner_red_sq;

        if (dist_sq_diff >= 0.0)
        {
          coll.emplace_back(7, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(7, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
      else
      {
        const double R_outer_sum = o_outer_R + p_R;
        const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

        if (post_dist_sq > R_outer_sum_sq) continue;  // no collision

        p[i].AddCollNumObj();
        const double post_dist = sqrt(post_dist_sq);
        const Vec3 u_n_vec = r_post_pos / post_dist;
        const double overlap = R_outer_sum - post_dist;

        const Vec3 r_pos = o_pos - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - R_outer_sum_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(7, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(7, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::DetectHollowSphereCollImage(const double Image_Coef,
                                            const int Reserve) noexcept
{
  std::vector<CollisionInfo> coll;
  coll.reserve(Reserve);

  std::vector<HollowSphereInfo>& o = obj->GetHollowSphere();
  const int hollow_sphere_size = obj->GetHollowSphereSize();

  for (int k = 0; k < hollow_sphere_size; ++k)
  {
    if (!(o[k].isCollOn())) continue;

    const Vec3 o_h_velo = o[k].GetHStepVelo();
    const Vec3 o_pos = o[k].GetPos();
    const Vec3 o_post_pos = o[k].GetPostPos();
    const double o_inner_R = o[k].GetInnerRadius();
    const double o_outer_R = o[k].GetOuterRadius();
    const double o_middle_R_sq =
      (o_inner_R + o_outer_R) * (o_inner_R + o_outer_R) * 0.25;

    std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
    const int ptcl_stay_size = ptcl->GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 r_post_pos = o_post_pos - p[i].GetPostPos();
      const double post_dist_sq = r_post_pos.LengthSq();
      const double p_R = p[i].GetRadius();

      if (post_dist_sq <= o_middle_R_sq)
      {
        const double R_inner_red = o_inner_R - p_R;
        const double R_inner_red_sq = R_inner_red * R_inner_red;

        const double post_dist = sqrt(post_dist_sq);
        const Vec3 u_n_vec = -r_post_pos / post_dist;

        // no collision
        if (R_inner_red_sq > post_dist_sq)
        {
          const double I_dist = o_inner_R - post_dist;
          const double I =
            GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec * I);
          continue;
        }

        p[i].AddCollNumObj();
        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const double overlap = post_dist - R_inner_red;
        const Vec3 r_pos = o_pos - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - R_inner_red_sq;

        if (dist_sq_diff >= 0.0)
        {
          coll.emplace_back(7, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(7, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
      else
      {
        const double R_outer_sum = o_outer_R + p_R;
        const double R_outer_sum_sq = R_outer_sum * R_outer_sum;

        const double post_dist = sqrt(post_dist_sq);
        const Vec3 u_n_vec = r_post_pos / post_dist;

        // no collision
        if (post_dist_sq > R_outer_sum_sq)
        {
          const double I_dist = post_dist - o_outer_R;
          const double I =
            GetImageForceLen(I_dist, Image_Coef, p[i].GetChgSq());
          p[i].AddImageForce(u_n_vec * I);
          continue;
        }

        p[i].AddCollNumObj();
        const double I = GetImageForceLen(p_R, Image_Coef, p[i].GetChgSq());
        p[i].AddImageForce(u_n_vec * I);

        const double overlap = R_outer_sum - post_dist;
        const Vec3 r_pos = o_pos - p[i].GetPos();
        const double dist_sq_diff = r_pos.LengthSq() - R_outer_sum_sq;

        if (dist_sq_diff <= 0.0)
        {
          coll.emplace_back(7, i, k, 0.0, overlap, u_n_vec, o_h_velo);
        }
        else
        {
          const Vec3 r_h_velo = o_h_velo - p[i].GetHStepVelo();
          const double time = GetCollTime(r_pos, r_h_velo, dist_sq_diff);
          coll.emplace_back(7, i, k, time, overlap, u_n_vec, o_h_velo);
        }
      }
    }
  }

  if (coll.empty()) return;

#pragma omp critical
  {
    coll_info.insert(coll_info.end(), std::make_move_iterator(coll.begin()),
                     std::make_move_iterator(coll.end()));
  }
}

void Collision::SortCollInfo() noexcept
{
  std::sort(coll_info.begin(), coll_info.end());
}

void Collision::CalcColl(const HardSphereInfo& C) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int coll_info_size = GetCollInfoSize();

  for (int i = 0; i < coll_info_size; ++i)
  {
    if (coll_info[i].obj_type != 0)
    {
      const auto [id_A, overlap, u_n_vec, h_velo_B] = coll_info[i].GetCollObj();
      p[id_A].AddPostPos(u_n_vec * -overlap);

      const double R_A = p[id_A].GetRadius();
      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      const Vec3 surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 pt_velo = h_velo_B - p[id_A].GetHStepVelo() - surf_velo_A;

      const double m_A = p[id_A].GetMass();
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_obj, m_A);
      const double Jt_len = GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_obj, m_A);
      const double fri_len = std::abs(Jn_len * C.coef_fri_obj);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-p[id_A].GetMassInv() * J);

      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      h_ang_velo_A -= (R_A * p[id_A].GetInertialMomentInv()) * J_ang_velo;
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
    }
    else
    {
      const auto [id_A, id_B, overlap, u_n_vec] = coll_info[i].GetCollPtcl();
      const double m_A = p[id_A].GetMass();
      const double m_B = p[id_B].GetMass();
      const double overlap_ratio = 1.0 / (m_A + m_B) * overlap;
      p[id_A].AddPostPos(-m_B * overlap_ratio * u_n_vec);
      p[id_B].AddPostPos(m_A * overlap_ratio * u_n_vec);

      const double R_A = p[id_A].GetRadius();
      const double R_B = p[id_B].GetRadius();
      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      Vec3 h_ang_velo_B = p[id_B].GetHStepAngVelo();
      const Vec3 h_surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 h_surf_velo_B = R_B * h_ang_velo_B.Cross(u_n_vec);
      const Vec3 r_h_velo = p[id_B].GetHStepVelo() - p[id_A].GetHStepVelo();
      const Vec3 pt_velo = r_h_velo - h_surf_velo_A - h_surf_velo_B;

      const double m_inv_A = p[id_A].GetMassInv();
      const double m_inv_B = p[id_B].GetMassInv();
      const double m_inv_AB_inv = 1.0 / (m_inv_A + m_inv_B);
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_ptcl, m_inv_AB_inv);
      const double Jt_len =
        GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_ptcl, m_inv_AB_inv);
      const double fri_len = std::abs(Jn_len * C.coef_fri_ptcl);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-m_inv_A * J);
      p[id_B].AddHStepVelo(m_inv_B * J);

      const double i_m_inv_A = p[id_A].GetInertialMomentInv();
      const double i_m_inv_B = p[id_B].GetInertialMomentInv();
      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      h_ang_velo_A -= (R_A * i_m_inv_A) * J_ang_velo;
      h_ang_velo_B -= (R_B * i_m_inv_B) * J_ang_velo;
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
      p[id_B].SetHStepAngVelo(h_ang_velo_B);
    }
  }
}

void Collision::CalcCollAdh(const HardSphereInfo& C) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int collision_info_num = static_cast<int>(coll_info.size());

  for (int i = 0; i < collision_info_num; ++i)
  {
    if (coll_info[i].obj_type != 0)
    {
      const auto [id_A, overlap, u_n_vec, h_velo_B] = coll_info[i].GetCollObj();
      p[id_A].AddPostPos(u_n_vec * -overlap);

      const double R_A = p[id_A].GetRadius();
      p[id_A].AddAdhForce(C.coef_adh * R_A * u_n_vec);

      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      const Vec3 surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 pt_velo = h_velo_B - p[id_A].GetHStepVelo() - surf_velo_A;

      const double m_A = p[id_A].GetMass();
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_obj, m_A);
      const double Jt_len = GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_obj, m_A);
      const double fri_len = std::abs(Jn_len * C.coef_fri_obj);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-p[id_A].GetMassInv() * J);

      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      h_ang_velo_A -= (R_A * p[id_A].GetInertialMomentInv()) * J_ang_velo;
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
    }
    else
    {
      const auto [id_A, id_B, overlap, u_n_vec] = coll_info[i].GetCollPtcl();
      const double m_A = p[id_A].GetMass();
      const double m_B = p[id_B].GetMass();
      const double overlap_ratio = 1.0 / (m_A + m_B) * overlap;
      p[id_A].AddPostPos(-m_B * overlap_ratio * u_n_vec);
      p[id_B].AddPostPos(m_A * overlap_ratio * u_n_vec);

      const double R_A = p[id_A].GetRadius();
      const double R_B = p[id_B].GetRadius();
      const Vec3 a_force = C.coef_adh * R_A * R_B / (R_A + R_B) * u_n_vec;
      p[id_A].AddAdhForce(a_force);
      p[id_B].AddAdhForce(-a_force);

      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      Vec3 h_ang_velo_B = p[id_B].GetHStepAngVelo();
      const Vec3 h_surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 h_surf_velo_B = R_B * h_ang_velo_B.Cross(u_n_vec);
      const Vec3 r_h_velo = p[id_B].GetHStepVelo() - p[id_A].GetHStepVelo();
      const Vec3 pt_velo = r_h_velo - h_surf_velo_A - h_surf_velo_B;

      const double m_inv_A = p[id_A].GetMassInv();
      const double m_inv_B = p[id_B].GetMassInv();
      const double m_inv_AB_inv = 1.0 / (m_inv_A + m_inv_B);
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_ptcl, m_inv_AB_inv);
      const double Jt_len =
        GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_ptcl, m_inv_AB_inv);
      const double fri_len = std::abs(Jn_len * C.coef_fri_ptcl);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-m_inv_A * J);
      p[id_B].AddHStepVelo(m_inv_B * J);

      const double i_m_inv_A = p[id_A].GetInertialMomentInv();
      const double i_m_inv_B = p[id_B].GetInertialMomentInv();
      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      h_ang_velo_A -= (R_A * i_m_inv_A) * J_ang_velo;
      h_ang_velo_B -= (R_B * i_m_inv_B) * J_ang_velo;
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
      p[id_B].SetHStepAngVelo(h_ang_velo_B);
    }
  }
}

void Collision::CalcCollRoll(const HardSphereInfo& C) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int collision_info_num = static_cast<int>(coll_info.size());

  for (int i = 0; i < collision_info_num; ++i)
  {
    if (coll_info[i].obj_type != 0)
    {
      const auto [id_A, overlap, u_n_vec, h_velo_B] = coll_info[i].GetCollObj();
      p[id_A].AddPostPos(u_n_vec * -overlap);

      const double R_A = p[id_A].GetRadius();
      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      const Vec3 surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 pt_velo = h_velo_B - p[id_A].GetHStepVelo() - surf_velo_A;

      const double m_A = p[id_A].GetMass();
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_obj, m_A);
      const double Jt_len = GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_obj, m_A);
      const double fri_len = std::abs(Jn_len * C.coef_fri_obj);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-p[id_A].GetMassInv() * J);

      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      const double i_m_inv_A = p[id_A].GetInertialMomentInv();
      h_ang_velo_A -= (R_A * i_m_inv_A) * J_ang_velo;

      const double roll_fri_coef_A = C.coef_roll_fri * R_A * i_m_inv_A;
      h_ang_velo_A -= GetRollFri(h_ang_velo_A, Jn_len, roll_fri_coef_A);
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
    }
    else
    {
      const auto [id_A, id_B, overlap, u_n_vec] = coll_info[i].GetCollPtcl();
      const double m_A = p[id_A].GetMass();
      const double m_B = p[id_B].GetMass();
      const double overlap_ratio = 1.0 / (m_A + m_B) * overlap;
      p[id_A].AddPostPos(-m_B * overlap_ratio * u_n_vec);
      p[id_B].AddPostPos(m_A * overlap_ratio * u_n_vec);

      const double R_A = p[id_A].GetRadius();
      const double R_B = p[id_B].GetRadius();
      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      Vec3 h_ang_velo_B = p[id_B].GetHStepAngVelo();
      const Vec3 h_surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 h_surf_velo_B = R_B * h_ang_velo_B.Cross(u_n_vec);
      const Vec3 r_h_velo = p[id_B].GetHStepVelo() - p[id_A].GetHStepVelo();
      const Vec3 pt_velo = r_h_velo - h_surf_velo_A - h_surf_velo_B;

      const double m_inv_A = p[id_A].GetMassInv();
      const double m_inv_B = p[id_B].GetMassInv();
      const double m_inv_AB_inv = 1.0 / (m_inv_A + m_inv_B);
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_ptcl, m_inv_AB_inv);
      const double Jt_len =
        GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_ptcl, m_inv_AB_inv);
      const double fri_len = std::abs(Jn_len * C.coef_fri_ptcl);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-m_inv_A * J);
      p[id_B].AddHStepVelo(m_inv_B * J);

      const double i_m_inv_A = p[id_A].GetInertialMomentInv();
      const double i_m_inv_B = p[id_B].GetInertialMomentInv();
      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      h_ang_velo_A -= (R_A * i_m_inv_A) * J_ang_velo;
      h_ang_velo_B -= (R_B * i_m_inv_B) * J_ang_velo;

      const double roll_fri_coef_A = C.coef_roll_fri * R_A * i_m_inv_A;
      const double roll_fri_coef_B = C.coef_roll_fri * R_B * i_m_inv_B;
      h_ang_velo_A -= GetRollFri(h_ang_velo_A, Jn_len, roll_fri_coef_A);
      h_ang_velo_B -= GetRollFri(h_ang_velo_B, Jn_len, roll_fri_coef_B);
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
      p[id_B].SetHStepAngVelo(h_ang_velo_B);
    }
  }
}

void Collision::CalcCollAdhRoll(const HardSphereInfo& C) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int collision_info_num = static_cast<int>(coll_info.size());

  for (int i = 0; i < collision_info_num; ++i)
  {
    if (coll_info[i].obj_type != 0)
    {
      const auto [id_A, overlap, u_n_vec, h_velo_B] = coll_info[i].GetCollObj();
      p[id_A].AddPostPos(u_n_vec * -overlap);

      const double R_A = p[id_A].GetRadius();
      p[id_A].AddAdhForce(C.coef_adh * R_A * u_n_vec);

      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      const Vec3 surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 pt_velo = h_velo_B - p[id_A].GetHStepVelo() - surf_velo_A;

      const double m_A = p[id_A].GetMass();
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_obj, m_A);
      const double Jt_len = GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_obj, m_A);
      const double fri_len = std::abs(Jn_len * C.coef_fri_obj);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-p[id_A].GetMassInv() * J);

      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      const double i_m_inv_A = p[id_A].GetInertialMomentInv();
      h_ang_velo_A -= (R_A * i_m_inv_A) * J_ang_velo;

      const double roll_fri_coef_A = C.coef_roll_fri * R_A * i_m_inv_A;
      h_ang_velo_A -= GetRollFri(h_ang_velo_A, Jn_len, roll_fri_coef_A);
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
    }
    else
    {
      const auto [id_A, id_B, overlap, u_n_vec] = coll_info[i].GetCollPtcl();
      const double m_A = p[id_A].GetMass();
      const double m_B = p[id_B].GetMass();
      const double overlap_ratio = 1.0 / (m_A + m_B) * overlap;
      p[id_A].AddPostPos(-m_B * overlap_ratio * u_n_vec);
      p[id_B].AddPostPos(m_A * overlap_ratio * u_n_vec);

      const double R_A = p[id_A].GetRadius();
      const double R_B = p[id_B].GetRadius();
      const Vec3 a_force = C.coef_adh * R_A * R_B / (R_A + R_B) * u_n_vec;
      p[id_A].AddAdhForce(a_force);
      p[id_B].AddAdhForce(-a_force);

      Vec3 h_ang_velo_A = p[id_A].GetHStepAngVelo();
      Vec3 h_ang_velo_B = p[id_B].GetHStepAngVelo();
      const Vec3 h_surf_velo_A = R_A * h_ang_velo_A.Cross(u_n_vec);
      const Vec3 h_surf_velo_B = R_B * h_ang_velo_B.Cross(u_n_vec);
      const Vec3 r_h_velo = p[id_B].GetHStepVelo() - p[id_A].GetHStepVelo();
      const Vec3 pt_velo = r_h_velo - h_surf_velo_A - h_surf_velo_B;

      const double m_inv_A = p[id_A].GetMassInv();
      const double m_inv_B = p[id_B].GetMassInv();
      const double m_inv_AB_inv = 1.0 / (m_inv_A + m_inv_B);
      const auto [Jn_len, u_t_vec] =
        GetJnLen(pt_velo, u_n_vec, C.coef_rest_n_ptcl, m_inv_AB_inv);
      const double Jt_len =
        GetJtLen(pt_velo, u_t_vec, C.coef_rest_t_ptcl, m_inv_AB_inv);
      const double fri_len = std::abs(Jn_len * C.coef_fri_ptcl);
      const Vec3 J = GetJ(fri_len, Jn_len, u_n_vec, Jt_len, u_t_vec);
      p[id_A].AddHStepVelo(-m_inv_A * J);
      p[id_B].AddHStepVelo(m_inv_B * J);

      const double i_m_inv_A = p[id_A].GetInertialMomentInv();
      const double i_m_inv_B = p[id_B].GetInertialMomentInv();
      const Vec3 J_ang_velo = u_n_vec.Cross(J);
      h_ang_velo_A -= (R_A * i_m_inv_A) * J_ang_velo;
      h_ang_velo_B -= (R_B * i_m_inv_B) * J_ang_velo;

      const double roll_fri_coef_A = C.coef_roll_fri * R_A * i_m_inv_A;
      const double roll_fri_coef_B = C.coef_roll_fri * R_B * i_m_inv_B;
      h_ang_velo_A -= GetRollFri(h_ang_velo_A, Jn_len, roll_fri_coef_A);
      h_ang_velo_B -= GetRollFri(h_ang_velo_B, Jn_len, roll_fri_coef_B);
      p[id_A].SetHStepAngVelo(h_ang_velo_A);
      p[id_B].SetHStepAngVelo(h_ang_velo_B);
    }
  }
}

}  // namespace dem
