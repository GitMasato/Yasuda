#include "Particle.h"

namespace dem
{
void Particle::SortPtclStay() noexcept
{
  std::sort(ptcl_stay.begin(), ptcl_stay.end());
}

void Particle::SortPtclOut() noexcept
{
  std::sort(ptcl_out.begin(), ptcl_out.end());
}

void Particle::LoadPtcl(const std::string& File_Name) noexcept
{
  // open ptcl data
  std::cout << "loading ptcl [" << File_Name << "] data..." << std::endl;
  std::ifstream fin(File_Name);
  if (!fin) data_io.ErrorInput(File_Name);

  std::string one_line;
  std::list<std::string> variable_list;
  std::list<std::string>::iterator iter;

  getline(fin, one_line);
  variable_list = data_io.Split(one_line, ",");
  iter = variable_list.begin();

  const int ptcl_size = data_io.StringToInt(*iter);
  ptcl_stay.resize(ptcl_size);
  getline(fin, one_line);  // ignore one line

  for (auto& v : ptcl_stay)
  {
    const int i = static_cast<int>(&v - &ptcl_stay[0]);
    v.SetID(i);
    getline(fin, one_line);
    variable_list = data_io.Split(one_line, ",");
    iter = variable_list.begin();
    Vec3 pos, velo, ang_velo;

    pos.x = data_io.StringToDouble(*iter);
    ++iter;
    pos.y = data_io.StringToDouble(*iter);
    ++iter;
    pos.z = data_io.StringToDouble(*iter);
    v.SetPostPos(pos);
    v.SetPos(pos);

    ++iter;
    velo.x = data_io.StringToDouble(*iter);
    ++iter;
    velo.y = data_io.StringToDouble(*iter);
    ++iter;
    velo.z = data_io.StringToDouble(*iter);
    v.SetHStepVelo(velo);
    v.SetVelo(velo);

    ++iter;
    ang_velo.x = data_io.StringToDouble(*iter);
    ++iter;
    ang_velo.y = data_io.StringToDouble(*iter);
    ++iter;
    ang_velo.z = data_io.StringToDouble(*iter);
    v.SetHStepAngVelo(ang_velo);
    v.SetAngVelo(ang_velo);

    ++iter;
    const double density = data_io.StringToDouble(*iter);
    ++iter;
    v.SetRadius(0.5 * data_io.StringToDouble(*iter));
    v.SetMass(4.0 / 3.0 * PI * v.GetRadiusCu() * density);
    v.SetMassInv(1.0 / v.GetMass());

    v.SetInertialMoment(2.0 / 5.0 * v.GetMass() * v.GetRadiusSq());
    v.SetInertialMomentInv(1.0 / v.GetInertialMoment());

    ++iter;
    v.SetChg(data_io.StringToDouble(*iter));
  }

  double min_radius = ptcl_stay[0].GetRadius();
  double max_radius = ptcl_stay[0].GetRadius();

  for (const auto& v : ptcl_stay)
  {
    const double radius = v.GetRadius();
    if (radius < min_radius) min_radius = radius;
    if (radius > max_radius) max_radius = radius;
  }

  SetMinDiameter(min_radius * 2.0);
  SetMaxDiameter(max_radius * 2.0);

  fin.close();
  std::cout << "particle number (" << GetPtclStaySize() << ")" << std::endl;
  std::cout << "finish loading ptcl [" << File_Name << "] data" << std::endl;
  std::cout << std::endl;
}

void Particle::CheckBoundaryX(const double DomainX) noexcept
{
  std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
  const double o_pos_x_1 = o[0].GetPosX();
  const double o_pos_x_2 = o_pos_x_1 + DomainX;
  const int ptcl_num = GetPtclStaySize();
  int counter = 0;

  for (int i = 0; i < ptcl_num - counter; i++)
  {
    const int t = i - counter;
    const double p_pos_x = ptcl_stay[t].GetPostPosX();

    if ((p_pos_x < o_pos_x_1) || (o_pos_x_2 <= p_pos_x))
    {
      ptcl_out.push_back(ptcl_stay[t]);
      ptcl_stay.erase(ptcl_stay.begin() + t);
      counter++;
    }
  }

  if (counter) is_ptcl_stay_size_changed = true;
}

void Particle::CheckBoundaryY(const double DomainY) noexcept
{
  std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
  const double o_pos_y_1 = o[0].GetPosY();
  const double o_pos_y_2 = o_pos_y_1 + DomainY;
  const int ptcl_num = GetPtclStaySize();
  int counter = 0;

  for (int i = 0; i < ptcl_num - counter; i++)
  {
    const int t = i - counter;
    const double p_pos_y = ptcl_stay[t].GetPostPosY();

    if ((p_pos_y < o_pos_y_1) || (o_pos_y_2 <= p_pos_y))
    {
      ptcl_out.push_back(ptcl_stay[t]);
      ptcl_stay.erase(ptcl_stay.begin() + t);
      counter++;
    }
  }

  if (counter) is_ptcl_stay_size_changed = true;
}

void Particle::CheckBoundaryZ(const double DomainZ) noexcept
{
  std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
  const double o_pos_z_1 = o[0].GetPosZ();
  const double o_pos_z_2 = o_pos_z_1 + DomainZ;
  const int ptcl_num = GetPtclStaySize();
  int counter = 0;

  for (int i = 0; i < ptcl_num - counter; i++)
  {
    const int t = i - counter;
    const double p_pos_z = ptcl_stay[t].GetPostPosZ();

    if ((p_pos_z < o_pos_z_1) || (o_pos_z_2 <= p_pos_z))
    {
      ptcl_out.push_back(ptcl_stay[t]);
      ptcl_stay.erase(ptcl_stay.begin() + t);
      counter++;
    }
  }

  if (counter) is_ptcl_stay_size_changed = true;
}

void Particle::CheckBoundaryXX(const double DomainX) noexcept
{
  std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
  const double o_pos_x_1 = o[0].GetPosX();
  const double o_pos_x_2 = o_pos_x_1 + DomainX;
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; i++)
  {
    double p_pos_x = ptcl_stay[i].GetPostPosX();

    if (p_pos_x < o_pos_x_1)
    {
      p_pos_x = fmod((p_pos_x - o_pos_x_1), DomainX) + o_pos_x_2;
    }
    else if (o_pos_x_2 <= p_pos_x)
    {
      p_pos_x = fmod((p_pos_x - o_pos_x_1), DomainX) + o_pos_x_1;
    }

    ptcl_stay[i].SetPostPosX(p_pos_x);
  }
}

void Particle::CheckBoundaryYY(const double DomainY) noexcept
{
  std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
  const double o_pos_y_1 = o[0].GetPosY();
  const double o_pos_y_2 = o_pos_y_1 + DomainY;
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; i++)
  {
    double p_pos_y = ptcl_stay[i].GetPostPosY();

    if (p_pos_y < o_pos_y_1)
    {
      p_pos_y = fmod((p_pos_y - o_pos_y_1), DomainY) + o_pos_y_2;
    }
    else if (o_pos_y_2 <= p_pos_y)
    {
      p_pos_y = fmod((p_pos_y - o_pos_y_1), DomainY) + o_pos_y_1;
    }

    ptcl_stay[i].SetPostPosY(p_pos_y);
  }
}

void Particle::CheckBoundaryZZ(const double DomainZ) noexcept
{
  std::vector<ObjectInfo>& o = obj->GetBoundaryBase();
  const double o_pos_z_1 = o[0].GetPosZ();
  const double o_pos_z_2 = o_pos_z_1 + DomainZ;
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; i++)
  {
    double p_pos_z = ptcl_stay[i].GetPostPosZ();

    if (p_pos_z < o_pos_z_1)
    {
      p_pos_z = fmod((p_pos_z - o_pos_z_1), DomainZ) + o_pos_z_2;
    }
    else if (o_pos_z_2 <= p_pos_z)
    {
      p_pos_z = fmod((p_pos_z - o_pos_z_1), DomainZ) + o_pos_z_1;
    }

    ptcl_stay[i].SetPostPosZ(p_pos_z);
  }
}

void Particle::ClearAdhForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; i++)
  {
    ptcl_stay[i].SetAdhForce(Vec3::Zero());
  }
}

void Particle::ClearImageForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; i++)
  {
    ptcl_stay[i].SetImageForce(Vec3::Zero());
  }
}

void Particle::ClearCoulombPtclForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; i++)
  {
    ptcl_stay[i].SetCoulombPtclForce(Vec3::Zero());
  }
}

void Particle::ClearCoulombExtForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].SetCoulombExtForce(Vec3::Zero());
  }
}

void Particle::ClearDieleForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].SetDieleForce(Vec3::Zero());
  }
}

void Particle::ClearDieleTorque() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].SetDieleTorque(Vec3::Zero());
  }
}

void Particle::ClearMagForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].SetMagForce(Vec3::Zero());
  }
}

void Particle::ClearMagTorque() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].SetMagTorque(Vec3::Zero());
  }
}

void Particle::StorePreExtForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].SetExtForce(ptcl_stay[i].GetPostExtForce());
    ptcl_stay[i].SetPostExtForce(Vec3::Zero());
  }
}

void Particle::AddGravityForce(const Vec3& Gravity) noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtForce(ptcl_stay[i].GetMass() * Gravity);
  }
}

void Particle::AddCoulombExtForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtForce(ptcl_stay[i].GetCoulombExtForce());
  }
}

void Particle::AddCoulombPtclForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtForce(ptcl_stay[i].GetCoulombPtclForce());
  }
}

void Particle::AddDieleForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtForce(ptcl_stay[i].GetDieleForce());
  }
}

void Particle::AddImageForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtForce(ptcl_stay[i].GetImageForce());
  }
}

void Particle::AddMagForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtForce(ptcl_stay[i].GetMagForce());
  }
}

void Particle::AddAdhForce() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtForce(ptcl_stay[i].GetAdhForce());
  }
}

void Particle::AddAirDrag(AirDragInfo& Info) noexcept
{
  const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const double air_coef = Info.coef_stokes * ptcl_stay[i].GetRadius();
    const Vec3 air_drag =
      -air_coef * (ptcl_stay[i].GetHStepVelo() - Info.air_velo);
    ptcl_stay[i].AddPostExtForce(air_drag);
  }
}

void Particle::AddAirDragCunningham(AirDragInfo& Info) noexcept
{
  if (Info.mean_free_path > 0.0)
  {
    const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const double R = ptcl_stay[i].GetRadius();
      const double Kn = Info.mean_free_path / (2.0 * R);
      const double Cc = 1.0 + Kn * (1.257 + 0.4 * exp(-1.1 / Kn));
      const double air_coef = Info.coef_stokes * R / Cc;
      const Vec3 air_drag =
        -air_coef * (ptcl_stay[i].GetHStepVelo() - Info.air_velo);
      ptcl_stay[i].AddPostExtForce(air_drag);
    }
  }
}

void Particle::AddAirDragReynolds(AirDragInfo& Info) noexcept
{
  if ((Info.coef_allen > 0.0) && (Info.coef_newton > 0.0) &&
      (Info.coef_stokes > 0.0))
  {
    const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
    for (int i = 0; i < ptcl_stay_size; ++i)
    {
      const Vec3 rel_velo = ptcl_stay[i].GetHStepVelo() - Info.air_velo;
      const double rel_velo_sq = rel_velo.LengthSq();
      const double R = ptcl_stay[i].GetRadius();
      const double reynolds_sq = Info.reynolds_coef_sq * rel_velo_sq * R * R;

      if (reynolds_sq < 4.0)
      {
        const double air_coef = Info.coef_stokes * R;
        const Vec3 air_drag = -air_coef * rel_velo;
        ptcl_stay[i].AddPostExtForce(air_drag);
      }
      else if (reynolds_sq < 250000.0)
      {
        const double reynolds_sqrt = pow(reynolds_sq, 0.25);
        const double air_coef = Info.coef_allen * R * R / reynolds_sqrt;
        const Vec3 velo_u_vec = rel_velo / sqrt(rel_velo_sq);
        const Vec3 air_drag = (-air_coef * rel_velo_sq) * velo_u_vec;
        ptcl_stay[i].AddPostExtForce(air_drag);
      }
      else
      {
        const double air_coef = Info.coef_newton * R * R;
        const Vec3 velo_u_vec = rel_velo / sqrt(rel_velo_sq);
        const Vec3 air_drag = (-air_coef * rel_velo_sq) * velo_u_vec;
        ptcl_stay[i].AddPostExtForce(air_drag);
      }
    }
  }
}

void Particle::StorePreExtTorque() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].SetExtTorque(ptcl_stay[i].GetPostExtTorque());
    ptcl_stay[i].SetPostExtTorque(Vec3::Zero());
  }
}

void Particle::AddDieleTorque() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtTorque(ptcl_stay[i].GetDieleTorque());
  }
}

void Particle::AddMagTorque() noexcept
{
  const int ptcl_num = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_num; ++i)
  {
    ptcl_stay[i].AddPostExtTorque(ptcl_stay[i].GetMagTorque());
  }
}

void Particle::TimeIntegHVeloPosVerlet(const double Time_Step) noexcept
{
  const double h_time = Time_Step * 0.5;
  const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const Vec3 accel = ptcl_stay[i].GetExtForce() * ptcl_stay[i].GetMassInv();
    const Vec3 h_step_velo = ptcl_stay[i].GetVelo() + accel * h_time;
    ptcl_stay[i].SetHStepVelo(h_step_velo);
    ptcl_stay[i].SetPostPos(ptcl_stay[i].GetPos() + h_step_velo * Time_Step);
  }
}

void Particle::TimeIntegVeloVerlet(const double Time_Step) noexcept
{
  const double h_time_step = Time_Step * 0.5;
  const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const Vec3 post_accel =
      ptcl_stay[i].GetPostExtForce() * ptcl_stay[i].GetMassInv();
    const Vec3 h_step_velo = ptcl_stay[i].GetHStepVelo();
    const Vec3 post_velo = h_step_velo + (post_accel * h_time_step);

    ptcl_stay[i].SetVelo(post_velo);
    ptcl_stay[i].SetPos(ptcl_stay[i].GetPostPos());
  }
}

void Particle::TimeIntegHAngVeloAngVerlet(const double Time_Step) noexcept
{
  const double h_time = Time_Step * 0.5;
  const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const Vec3 ang_accel =
      ptcl_stay[i].GetExtTorque() * ptcl_stay[i].GetInertialMomentInv();
    const Vec3 h_step_ang_velo = ptcl_stay[i].GetAngVelo() + ang_accel * h_time;
    ptcl_stay[i].SetHStepAngVelo(h_step_ang_velo);
    ptcl_stay[i].AddAng(h_step_ang_velo * Time_Step);
  }
}

void Particle::TimeIntegAngVeloVerlet(const double Time_Step) noexcept
{
  const double h_time_step = Time_Step * 0.5;
  const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const Vec3 post_ext_torque = ptcl_stay[i].GetPostExtTorque();
    const Vec3 post_ang_accel =
      post_ext_torque * ptcl_stay[i].GetInertialMomentInv();
    const Vec3 h_step_ang_velo = ptcl_stay[i].GetHStepAngVelo();
    const Vec3 post_ang_velo = h_step_ang_velo + (post_ang_accel * h_time_step);
    ptcl_stay[i].SetAngVelo(post_ang_velo);
  }
}

std::array<double, 2> Particle::GetKineticEnergy() noexcept
{
#pragma omp single
  {
    trans_energy = 0.0;
    rot_energy = 0.0;
  }

  double t_energy = 0.0;
  double r_energy = 0.0;
  const int ptcl_stay_size = GetPtclStaySize();

#pragma omp for schedule(dynamic)
  for (int i = 0; i < ptcl_stay_size; ++i)
  {
    const Vec3 velo = ptcl_stay[i].GetVelo();
    t_energy += 0.5 * ptcl_stay[i].GetMass() * velo.LengthSq();
    const Vec3 ang_velo = ptcl_stay[i].GetAngVelo();
    r_energy += 0.5 * ptcl_stay[i].GetInertialMoment() * ang_velo.LengthSq();
  }

#pragma omp atomic
  trans_energy += t_energy;
#pragma omp atomic
  rot_energy += r_energy;
#pragma omp barrier

  return {trans_energy, rot_energy};
}

}  // namespace dem
