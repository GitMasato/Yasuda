#pragma once
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Clock.h"
#include "DataIO.h"
#include "Field.h"
#include "Particle.h"

namespace dem
{
class Output
{
private:
  Field* ele_field;
  Field* mag_field;
  Particle* ptcl;
  Clock clock;
  DataIO data_io;
  constexpr static double PI = 3.14159;

public:
  Output() noexcept
    : ele_field(nullptr), mag_field(nullptr), ptcl(nullptr)
  {
  }

  void AssociatePtcl(Particle* Ptcl) noexcept { ptcl = Ptcl; }

  void AssociateEleField(Field* Ele_Field) noexcept { ele_field = Ele_Field; }

  void AssociateMagField(Field* Mag_Field) noexcept { mag_field = Mag_Field; }

  void CreateTimeFile(const int Total_Time_Step,
                      const std::array<double, 2> Energy,
                      const std::string& File) const noexcept;

  void WriteTimeFile(const int Time_Step, const int Total_Time_Step,
                     const std::array<double, 2> Energy,
                     const std::string& File) const noexcept;

  void WriteTimeFileEnd(const std::string& File) const noexcept;

  void CreatePtclStayFile(const std::string& File) const noexcept;

  void CreatePtclOutFile(const std::string& File) const noexcept;

  void CreatePtclFile(const double Elapsed_Time, const std::string& File,
                      std::vector<std::string> Keys) const noexcept;

  void WriteID(const char Key[3], std::ofstream& Fout,
               std::vector<ParticleInfo>& P) const noexcept;

  void WriteCollNumPtcl(const char Key[3], std::ofstream& Fout,
                        std::vector<ParticleInfo>& P) const noexcept;

  void WriteCollNumObj(const char Key[3], std::ofstream& Fout,
                       std::vector<ParticleInfo>& P) const noexcept;

  void WriteRadius(const char Key[3], std::ofstream& Fout,
                   std::vector<ParticleInfo>& P) const noexcept;

  void WriteMass(const char Key[3], std::ofstream& Fout,
                 std::vector<ParticleInfo>& P) const noexcept;

  void WriteUnitCharge(const char Key[3], std::ofstream& Fout,
                       std::vector<ParticleInfo>& P) const noexcept;

  void WriteCharge(const char Key[3], std::ofstream& Fout,
                   std::vector<ParticleInfo>& P) const noexcept;

  void WriteEleDipoleEnergy(const char Key[3], std::ofstream& Fout) const
    noexcept;

  void WriteMagDipoleEnergy(const char Key[3], std::ofstream& Fout) const
    noexcept;

  void WritePos(const char Key[3], std::ofstream& Fout,
                std::vector<ParticleInfo>& P) const noexcept;

  void WriteVelo(const char Key[3], std::ofstream& Fout,
                 std::vector<ParticleInfo>& P) const noexcept;

  void WriteAngVelo(const char Key[3], std::ofstream& Fout,
                    std::vector<ParticleInfo>& P) const noexcept;

  void WriteAng(const char Key[3], std::ofstream& Fout,
                std::vector<ParticleInfo>& P) const noexcept;

  void WriteForce(const char Key[3], std::ofstream& Fout,
                  std::vector<ParticleInfo>& P) const noexcept;

  void WriteTorque(const char Key[3], std::ofstream& Fout,
                   std::vector<ParticleInfo>& P) const noexcept;

  void WriteCoulombExtForce(const char Key[3], std::ofstream& Fout,
                            std::vector<ParticleInfo>& P) const noexcept;

  void WriteCoulombPtclForce(const char Key[3], std::ofstream& Fout,
                             std::vector<ParticleInfo>& P) const noexcept;

  void WriteImageForce(const char Key[3], std::ofstream& Fout,
                       std::vector<ParticleInfo>& P) const noexcept;

  void WriteDieleForce(const char Key[3], std::ofstream& Fout,
                       std::vector<ParticleInfo>& P) const noexcept;

  void WriteMagForce(const char Key[3], std::ofstream& Fout,
                     std::vector<ParticleInfo>& P) const noexcept;

  void WriteAdhForce(const char Key[3], std::ofstream& Fout,
                     std::vector<ParticleInfo>& P) const noexcept;

  void WriteDieleTorque(const char Key[3], std::ofstream& Fout,
                        std::vector<ParticleInfo>& P) const noexcept;

  void WriteMagTorque(const char Key[3], std::ofstream& Fout,
                      std::vector<ParticleInfo>& P) const noexcept;
};

}  // namespace dem
