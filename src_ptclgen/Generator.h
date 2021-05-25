#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <new>
#include <numeric>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "DataIO.h"
#include "ParticleInfo.h"
#include "Vector.h"

namespace dem {

class Generator {
private:
  DataIO data_io;
  std::vector<ParticleInfo> ptcl_data;
  std::string gen_file;
  std::string ptcl_file;
  std::vector<double> pos_params;
  std::vector<double> velo_params;
  std::vector<double> ang_velo_params;
  std::vector<double> density_params;
  std::vector<double> diameter_params;
  std::vector<double> charge_params;
  int ptcl_size;
  int pos_setting_type, velo_setting_type, ang_velo_setting_type;
  int density_setting_type, diameter_setting_type, charge_setting_type;

public:
  Generator() = default;
  ~Generator(){};

  void read_args(const int &argc, char **argv);
  void read_gen() noexcept;
  void read_ptcl() noexcept;
  void set();
  void set_pos();
  void set_velo();
  void set_ang_velo();
  void set_density();
  void set_diameter();
  void set_charge();
  void generate();

  template <typename T>
  void set_variable(const std::vector<std::string> &Value, const std::string &Key,
                    T &Parameter) noexcept {
    if (Value[0] == Key) {
      if (Value.size() == 2) {
        data_io.TransformString(Value[1], Parameter);
      } else {
        data_io.ErrorInput(Key);
      }
    }
  }

  template <typename T>
  void set_parameter(const std::vector<std::string> &Value, const std::string &Key,
                     int &Type, T &Parameter) noexcept {
    if (Value[0] == Key) {
      const std::size_t values = Value.size();

      if (values >= 3) {
        data_io.TransformString(Value[1], Type);

        if (Type) {
          Parameter.resize(values - 2);

          for (std::size_t i = 2; i < values; ++i) {
            if (Type) data_io.TransformString(Value[i], Parameter[i - 2]);
          }
        }
      } else {
        data_io.ErrorInput(Key);
      }
    }
  }
};

}  // namespace dem
