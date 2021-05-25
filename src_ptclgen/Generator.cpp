#include "Generator.h"

namespace dem {

void Generator::read_args(const int &argc, char **argv) {

  if (argc <= 1) {
    std::cout << "no input file given..." << std::endl;
    std::exit(1);
  }

  std::filesystem::path input = argv[1];
  gen_file =
    ((input.is_absolute()) ? input : std::filesystem::current_path() / input).string();
  read_gen();

  if (argc == 2) {

    if (!ptcl_size) {
      std::cout << "ptcl number is 0!" << std::endl;
      std::exit(1);
    }

    ptcl_data.resize(ptcl_size);

  } else if (argc == 3) {

    input = argv[2];
    ptcl_file =
      ((input.is_absolute()) ? input : std::filesystem::current_path() / input)
        .string();
    read_ptcl();

  } else {
    std::cout << "too much arguments!" << std::endl;
    std::exit(1);
  }
}

void Generator::read_gen() noexcept {

  if (!std::filesystem::exists(gen_file)) {
    std::cout << "'" << gen_file << "' does not exist!" << std::endl;
    std::exit(1);
  }

  if (gen_file.find(".ptclgen") == std::string::npos) {
    std::cout << "'" << gen_file << "' is not '.ptclgen' file!" << std::endl;
    std::exit(1);
  }

  std::cout << "read '" << gen_file << "' ..." << std::endl;
  std::ifstream fin(gen_file);
  if (!fin) { data_io.ErrorInput(gen_file); }

  std::string line;

  while (!fin.eof() && getline(fin, line)) {

    std::vector<std::string> key;
    data_io.SplitLine(line, key);
    if (key.empty()) continue;

    set_parameter(key, "pos", pos_setting_type, pos_params);
    set_parameter(key, "velo", velo_setting_type, velo_params);
    set_parameter(key, "ang_velo", ang_velo_setting_type, ang_velo_params);
    set_parameter(key, "density", density_setting_type, density_params);
    set_parameter(key, "diameter", diameter_setting_type, diameter_params);
    set_parameter(key, "charge", charge_setting_type, charge_params);
    set_variable(key, "particle_number", ptcl_size);
  }
}

void Generator::read_ptcl() noexcept {

  if (!std::filesystem::exists(ptcl_file)) {
    std::cout << "'" << ptcl_file << "' does not exist!" << std::endl;
    std::exit(1);
  }

  if (ptcl_file.find("_ptcl.csv") == std::string::npos) {
    std::cout << "'" << gen_file << "' is not '_ptcl.csv' file!" << std::endl;
    std::exit(1);
  }

  std::cout << "read '" << ptcl_file << "' ..." << std::endl;
  std::ifstream fin(ptcl_file);
  if (!fin) { data_io.ErrorInput(ptcl_file); }

  std::string one_line;
  std::list<std::string> variable_list;
  std::list<std::string>::iterator iter;

  getline(fin, one_line);
  variable_list = data_io.Split(one_line, ",");
  iter = variable_list.begin();
  ptcl_size = data_io.StringToInt(*iter);
  ptcl_data.clear();
  ptcl_data.resize(ptcl_size);
  getline(fin, one_line);

  for (auto &p : ptcl_data) {

    getline(fin, one_line);
    variable_list = data_io.Split(one_line, ",");
    iter = variable_list.begin();

    p.pos.x = data_io.StringToDouble(*iter);
    ++iter;
    p.pos.y = data_io.StringToDouble(*iter);
    ++iter;
    p.pos.z = data_io.StringToDouble(*iter);
    ++iter;
    p.velo.x = data_io.StringToDouble(*iter);
    ++iter;
    p.velo.y = data_io.StringToDouble(*iter);
    ++iter;
    p.velo.z = data_io.StringToDouble(*iter);
    ++iter;
    p.ang_velo.x = data_io.StringToDouble(*iter);
    ++iter;
    p.ang_velo.y = data_io.StringToDouble(*iter);
    ++iter;
    p.ang_velo.z = data_io.StringToDouble(*iter);
    ++iter;
    p.density = data_io.StringToDouble(*iter);
    ++iter;
    p.diameter = data_io.StringToDouble(*iter);
    ++iter;
    p.charge = data_io.StringToDouble(*iter);
  }
}

void Generator::set() {

  if (density_setting_type) set_density();
  if (diameter_setting_type) set_diameter();
  if (charge_setting_type) set_charge();

  if (velo_setting_type) set_velo();
  if (ang_velo_setting_type) set_ang_velo();
  if (pos_setting_type) set_pos();
}

void Generator::set_pos() {

  for (auto &p : ptcl_data) { p.pos.Set(0.0, 0.0, 0.0); }
  std::random_device rnddev;
  std::mt19937 mt(rnddev());

  switch (pos_setting_type) {

    case 1: {  // args: x_min x_max y_min y_max z_min z_max

      std::uniform_real_distribution<double> rnd_x(pos_params[0], pos_params[1]);
      std::uniform_real_distribution<double> rnd_y(pos_params[2], pos_params[3]);
      std::uniform_real_distribution<double> rnd_z(pos_params[4], pos_params[5]);

      for (int i = 0; i < ptcl_size; i++) {

      point_1:

        const Vec3 gen(rnd_x(mt), rnd_y(mt), rnd_z(mt));

        for (int j = i - 1; j >= 0; j--) {
          const double r_sum = (ptcl_data[i].diameter + ptcl_data[j].diameter) * 0.5;
          if ((ptcl_data[j].pos - gen).LengthSq() < r_sum * r_sum) goto point_1;
        }

        ptcl_data[i].pos = gen;
        std::cout << "remaining " << ptcl_size - i << '\r';
      }

      break;
    }

    case 2: {  // args: x_center y_center z_center radius

      const Vec3 domain_center(pos_params[0], pos_params[1], pos_params[2]);
      const Vec3 pos_min = domain_center - pos_params[3];
      const Vec3 pos_max = domain_center + pos_params[3];
      std::uniform_real_distribution<double> rnd_x(pos_min.x, pos_max.x);
      std::uniform_real_distribution<double> rnd_y(pos_min.y, pos_max.y);
      std::uniform_real_distribution<double> rnd_z(pos_min.z, pos_max.z);
      const double domain_r_sq = pos_params[3] * pos_params[3];

      for (int i = 0; i < ptcl_size; i++) {

      point_3:

        const Vec3 gen(rnd_x(mt), rnd_y(mt), rnd_z(mt));
        if ((gen - domain_center).LengthSq() > domain_r_sq) goto point_3;

        for (int j = i - 1; j >= 0; j--) {
          const double r_sum = (ptcl_data[i].diameter + ptcl_data[j].diameter) * 0.5;
          if ((ptcl_data[j].pos - gen).LengthSq() < r_sum * r_sum) goto point_3;
        }

        ptcl_data[i].pos = gen;
        std::cout << "remaining " << ptcl_size - i << '\r';
      }

      break;
    }

    default: break;
  }
}

void Generator::set_velo() {

  for (auto &p : ptcl_data) { p.velo.Set(0.0, 0.0, 0.0); }
  std::random_device rnddev;
  std::mt19937 mt(rnddev());

  switch (velo_setting_type) {

    case 1: {  // args: x_min x_max y_min y_max z_min z_max

      std::uniform_real_distribution<double> rnd_x(velo_params[0], velo_params[1]);
      std::uniform_real_distribution<double> rnd_y(velo_params[2], velo_params[3]);
      std::uniform_real_distribution<double> rnd_z(velo_params[4], velo_params[5]);
      for (auto &p : ptcl_data) { p.velo.Set(rnd_x(mt), rnd_y(mt), rnd_z(mt)); }
      break;
    }

    case 2: {  // args: x_ave x_dev y_ave y_dev z_ave z_dev

      std::normal_distribution<double> rnd_x(velo_params[0], velo_params[1]);
      std::normal_distribution<double> rnd_y(velo_params[2], velo_params[3]);
      std::normal_distribution<double> rnd_z(velo_params[4], velo_params[5]);
      for (auto &p : ptcl_data) { p.velo.Set(rnd_x(mt), rnd_y(mt), rnd_z(mt)); }
      break;
    }

    default: break;
  }
}

void Generator::set_ang_velo() {

  for (auto &p : ptcl_data) { p.ang_velo.Set(0.0, 0.0, 0.0); }
  std::random_device rnddev;
  std::mt19937 mt(rnddev());

  switch (ang_velo_setting_type) {

    case 1: {  // args: x_min x_max y_min y_max z_min z_max

      std::uniform_real_distribution<double> rnd_x(ang_velo_params[0],
                                                   ang_velo_params[1]);
      std::uniform_real_distribution<double> rnd_y(ang_velo_params[2],
                                                   ang_velo_params[3]);
      std::uniform_real_distribution<double> rnd_z(ang_velo_params[4],
                                                   ang_velo_params[5]);
      for (auto &p : ptcl_data) { p.ang_velo.Set(rnd_x(mt), rnd_y(mt), rnd_z(mt)); }
      break;
    }

    case 2: {  // args: x_ave x_dev y_ave y_dev z_ave z_dev

      std::normal_distribution<double> rnd_x(ang_velo_params[0], ang_velo_params[1]);
      std::normal_distribution<double> rnd_y(ang_velo_params[2], ang_velo_params[3]);
      std::normal_distribution<double> rnd_z(ang_velo_params[4], ang_velo_params[5]);
      for (auto &p : ptcl_data) { p.ang_velo.Set(rnd_x(mt), rnd_y(mt), rnd_z(mt)); }
      break;
    }

    default: break;
  }
}

void Generator::set_density() {

  for (auto &p : ptcl_data) { p.density = 0.0; }
  std::random_device rnddev;
  std::mt19937 mt(rnddev());

  switch (density_setting_type) {

    case 1: {  // args: min max

      std::uniform_real_distribution<double> rnd(density_params[0], density_params[1]);
      for (auto &p : ptcl_data) { p.density = rnd(mt); }
      break;
    }

    case 2: {  // args: ave dev

      std::normal_distribution<double> rnd(density_params[0], density_params[1]);
      double tmp;

      for (auto &p : ptcl_data) {
        do { tmp = rnd(mt); } while (tmp < 0);
        p.density = tmp;
      }

      break;
    }

    default: break;
  }
}

void Generator::set_diameter() {

  for (auto &p : ptcl_data) { p.diameter = 0.0; }
  std::random_device rnddev;
  std::mt19937 mt(rnddev());

  switch (diameter_setting_type) {

    case 1: {  // args: min max

      std::uniform_real_distribution<double> rnd(diameter_params[0],
                                                 diameter_params[1]);
      for (auto &p : ptcl_data) { p.diameter = rnd(mt); }
      break;
    }
    case 2: {  // args: ave dev

      std::normal_distribution<double> rnd(diameter_params[0], diameter_params[1]);
      double tmp;

      for (auto &p : ptcl_data) {
        do { tmp = rnd(mt); } while (tmp < 0);
        p.diameter = tmp;
      }

      break;
    }

    default: break;
  }

  std::sort(ptcl_data.begin(), ptcl_data.end(),
            [](const ParticleInfo &a, const ParticleInfo &b) {
              return a.diameter < b.diameter;
            });
}

void Generator::set_charge() {

  for (auto &p : ptcl_data) { p.charge = 0.0; }
  std::random_device rnddev;
  std::mt19937 mt(rnddev());

  switch (charge_setting_type) {

    case 1: {  // args: min max

      std::uniform_real_distribution<double> rnd(charge_params[0], charge_params[1]);
      for (auto &p : ptcl_data) { p.charge = rnd(mt); }
      break;
    }

    case 2: {  // args: ave dev

      std::normal_distribution<double> rnd(charge_params[0], charge_params[1]);
      for (auto &p : ptcl_data) { p.charge = rnd(mt); }
      break;
    }

    default: break;
  }
}

void Generator::generate() {

  if (ptcl_data.empty()) { std::cout << "no ptcl data is generated..." << std::endl; }

  namespace fs = std::filesystem;
  fs::path p_cr = fs::current_path();
  fs::path p_gen = gen_file;
  std::string gen_name = p_gen.filename().string();
  gen_name.erase(gen_name.end() - 8, gen_name.end());
  fs::path p_output = p_cr;
  p_output.append(gen_name + "_ptcl.csv");

  std::ofstream fout;
  fout.open(p_output.string(), std::ios_base::out | std::ios_base::trunc);
  std::cout << "generate '" << p_output.string() + "'..." << std::endl;

  fout << ptcl_data.size() << ",number" << std::endl;
  fout << "x (m),y (m),z (m),v_x (m/s),v_y (m/s),v_z (m/s),a_v_x (rad/s),a_v_y "
          "(rad/s),a_v_z (rad/s),density (kg/m3),diameter (m),charge (C)"
       << std::endl;

  for (auto &p : ptcl_data) {
    fout << p.pos.x << "," << p.pos.y << "," << p.pos.z << "," << p.velo.x << ","
         << p.velo.y << "," << p.velo.z << "," << p.ang_velo.x << "," << p.ang_velo.y
         << "," << p.ang_velo.z << "," << p.density << "," << p.diameter << ","
         << p.charge << std::endl;
  }

  fout.close();
}
}  // namespace dem
