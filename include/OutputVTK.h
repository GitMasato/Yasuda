#pragma once
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "DataIO.h"
#include "Field.h"
#include "Object.h"
#include "Particle.h"

namespace dem
{
class OutputVTK
{
private:
  Field* ele_field;
  Field* mag_field;
  Object* obj;
  Particle* ptcl;
  DataIO data_io;
  constexpr static double PI = 3.14159;

public:
  OutputVTK() noexcept
    : ele_field(nullptr), mag_field(nullptr), obj(nullptr), ptcl(nullptr)
  {
  }

  void AssociateObj(Object* Obj) noexcept { obj = Obj; }

  void AssociatePtcl(Particle* Ptcl) noexcept { ptcl = Ptcl; }

  void AssociateEleField(Field* Ele_Field) noexcept { ele_field = Ele_Field; }

  void AssociateMagField(Field* Mag_Field) noexcept { mag_field = Mag_Field; }

  void CreateFile(const double Elapsed_time, const std::string& File_Name,
                  std::vector<std::string> Keys) const noexcept;

  void WritePos(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteCellType(std::ofstream& Fout) const noexcept;

  void WriteID(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteCollNumPtcl(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteCollNumObj(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteRadius(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteMass(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteUnitCharge(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteCharge(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteEleDipoleEnergy(std::ofstream& Fout) const noexcept;

  void WriteMagDipoleEnergy(std::ofstream& Fout) const noexcept;

  void WriteVelo(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteAngVelo(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteAng(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteForce(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteTorque(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteCoulombExtForce(std::ofstream& Fout,
                            std::vector<ParticleInfo>& P) const noexcept;

  void WriteCoulombPtclForce(std::ofstream& Fout,
                             std::vector<ParticleInfo>& P) const noexcept;

  void WriteImageForce(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteDieleForce(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteMagForce(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteAdhForce(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteDieleTorque(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;

  void WriteMagTorque(std::ofstream& Fout, std::vector<ParticleInfo>& P) const
    noexcept;
};

}  // namespace dem
