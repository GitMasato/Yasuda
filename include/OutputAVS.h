#pragma once
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "DataIO.h"
#include "Object.h"
#include "Particle.h"

namespace dem
{
class OutputAVS
{
private:
  Object* obj;
  Particle* ptcl;
  DataIO data_io;
  int counter;
  constexpr static double PI = 3.14159;

public:
  OutputAVS() noexcept : obj(nullptr), ptcl(nullptr), counter(0) {}

  void AssociatePtcl(Particle* Ptcl) noexcept { ptcl = Ptcl; }

  void AssociateObj(Object* Obj) noexcept { obj = Obj; }

  void CreateAVS(const int Step, const Vec3& Gravity, const Vec3& Domain,
                 const std::string& Avs_File, const std::string& Key) noexcept;

  void WriteAVS(const double Time, const Vec3& Gravity, const Vec3& Domain,
                const std::string& Avs_File, const std::string& Key) noexcept;

  void WriteHeader(const double Time, std::ofstream& Fout) noexcept;

  void WritePtcl(std::ofstream& Fout) const noexcept;

  void WritePtclCharge(std::ofstream& Fout) const noexcept;

  void WritePtclForce(std::ofstream& Fout, const Vec3& Gravity) const noexcept;

  void WriteObj(std::ofstream& Fout, const Vec3& Domain) const noexcept;

  void WriteFloor(std::ofstream& Fout, const Vec3& Domain) const noexcept;

  void WriteBox(std::ofstream& Fout) const noexcept;

  void WriteHollowBox(std::ofstream& Fout) const noexcept;

  void WriteCylinder(std::ofstream& Fout) const noexcept;

  void WriteHollowCylinder(std::ofstream& Fout) const noexcept;

  void WriteSphere(std::ofstream& Fout) const noexcept;

  void WriteHollowSphere(std::ofstream& Fout) const noexcept;
};

}  // namespace dem
