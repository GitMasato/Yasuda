#include "LongRange.h"

namespace dem
{
void LongRange::CalcEleDipoleInteract(const double Dipole_Coef,
                                      const double Field_Coef) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_size = ptcl->GetPtclStaySize();

  std::vector<FieldInfoAtParticle>& f = field->GetFieldInfoAtPtcl();
  std::vector<LongRangePtclInfo>& lr = lr_ptcl_list->GetLRPtclInfo();

  do
  {
#pragma omp single
    {
      pre_sum_energy = sum_energy;
      sum_energy = 0;
    }

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_size; ++i)
    {
      Vec3 field_from_ptcl(0.0, 0.0, 0.0);

      for (const auto& v : lr[i].lr_ptcl)
      {
        const Vec3 rel_pos = v.rel_pos;
        const Vec3 vs_moment = f[v.id].GetMoment();
        const double d = vs_moment.Dot(rel_pos);
        const double coef_1 = 3.0 * Field_Coef * v.dist_inv_5;
        const double coef_2 = Field_Coef * v.dist_inv_3;
        field_from_ptcl += (coef_1 * d * rel_pos) - (coef_2 * vs_moment);
      }

      f[i].SetFieldPtcl(field_from_ptcl);
    }

    double dipole_energy = 0.0;

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_size; ++i)
    {
      const Vec3 field_total = f[i].GetFieldTotal();
      const Vec3 moment = p[i].GetRadiusCu() * Dipole_Coef * field_total;
      f[i].SetMoment(moment);
      dipole_energy += field_total.Dot(moment);
    }

#pragma omp critical
    sum_energy += dipole_energy;
#pragma omp barrier

  } while (fabs(sum_energy - pre_sum_energy) > EPS_E * ptcl_size);
}

void LongRange::CalcMagDipoleInteract(const double Dipole_Coef,
                                      const double Field_Coef) noexcept
{
  std::vector<ParticleInfo>& p = ptcl->GetPtclStay();
  const int ptcl_size = ptcl->GetPtclStaySize();

  std::vector<FieldInfoAtParticle>& f = field->GetFieldInfoAtPtcl();
  std::vector<LongRangePtclInfo>& lr = lr_ptcl_list->GetLRPtclInfo();

  do
  {
#pragma omp single
    {
      pre_sum_energy = sum_energy;
      sum_energy = 0;
    }

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_size; ++i)
    {
      Vec3 field_from_ptcl(0.0, 0.0, 0.0);

      for (const auto& v : lr[i].lr_ptcl)
      {
        const Vec3 rel_pos = v.rel_pos;
        const Vec3 vs_moment = f[v.id].GetMoment();
        const double d = vs_moment.Dot(rel_pos);
        const double coef_1 = 3.0 * Field_Coef * v.dist_inv_5;
        const double coef_2 = Field_Coef * v.dist_inv_3;
        field_from_ptcl += (coef_1 * d * rel_pos) - (coef_2 * vs_moment);
      }

      f[i].SetFieldPtcl(field_from_ptcl);
    }

    double dipole_energy = 0.0;

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_size; ++i)
    {
      const Vec3 total_field = f[i].GetFieldTotal();
      const Vec3 moment = p[i].GetRadiusCu() * Dipole_Coef * total_field;
      f[i].SetMoment(moment);
      dipole_energy += total_field.Dot(moment);
    }

#pragma omp critical
    sum_energy += dipole_energy;
#pragma omp barrier

  } while (fabs(sum_energy - pre_sum_energy) > EPS_B * ptcl_size);
}

void LongRange::SetFieldFromOtherDipole(const double Field_Coef) noexcept
{
  const int ptcl_stay_size = ptcl->GetPtclStaySize();
  std::vector<FieldInfoAtParticle>& f = field->GetFieldInfoAtPtcl();
  std::vector<LongRangePtclInfo>& lr = lr_ptcl_list->GetLRPtclInfo();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    std::array<Vec3, 3> fd{Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};

    for (const auto& v : lr[i].lr_ptcl)
    {
      const Vec3 rel_pos = v.rel_pos;
      const Vec3 vs_moment = f[v.id].GetMoment();
      const Vec3 Px = vs_moment.x * rel_pos;
      const Vec3 Py = vs_moment.y * rel_pos;
      const Vec3 Pz = vs_moment.z * rel_pos;
      const double P = Px.x + Py.y + Pz.z;
      const double C = 3.0 * Field_Coef * v.dist_inv_5;
      const double D = -15.0 * Field_Coef * P * v.dist_inv_7;

      const double xx = C * (2.0 * Px.x + P) + D * rel_pos.x * rel_pos.x;
      const double yy = C * (2.0 * Py.y + P) + D * rel_pos.y * rel_pos.y;
      const double zz = C * (2.0 * Pz.z + P) + D * rel_pos.z * rel_pos.z;
      const double xy = C * (Py.x + Px.y) + D * rel_pos.x * rel_pos.y;
      const double xz = C * (Pz.x + Px.z) + D * rel_pos.x * rel_pos.z;
      const double yz = C * (Pz.y + Py.z) + D * rel_pos.y * rel_pos.z;

      const Vec3 fdx(xx, xy, xz);
      const Vec3 fdy(xy, yy, yz);
      const Vec3 fdz(xz, yz, zz);
      fd[0] += fdx;
      fd[1] += fdy;
      fd[2] += fdz;
    }

    f[i].AddFieldDeriv(fd);
    f[i].SetFieldAsTotalField();
  }
}

}  // namespace dem
