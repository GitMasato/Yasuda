#pragma once

#include "Vector.h"

namespace dem
{

  struct ParticleInfo
  {
    Vec3 pos;
    Vec3 velo;
    Vec3 ang_velo;
    double density, diameter, charge;

    ParticleInfo() noexcept
        : pos(0.0, 0.0, 0.0), velo(0.0, 0.0, 0.0), ang_velo(0.0, 0.0, 0.0), density(0.0), diameter(0.0), charge(0.0) {}
  };
}
