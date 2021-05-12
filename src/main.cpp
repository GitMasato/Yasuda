#include <clipp/clipp.h>

#include "DEM.h"

int main(int argc, char** argv)
{
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  enum class mode
  {
    dem,
    ele,
    mag,
    help
  };
  mode selected = mode::help;

  auto dem_mode =
    "dem:" % (clipp::command("dem").set(selected, mode::dem),
              clipp::option("--replace") % "replace existing files",
              clipp::option("-f", "--force") % "don't ask for confirmation");

  auto ele_mode =
    "ele:" % (clipp::command("ele").set(selected, mode::ele),
              (clipp::command("date") | clipp::command("content")),
              clipp::option("-b", "--binary") % "compare files byte by byte",
              clipp::option("-q", "--quick") % "use heuristics for faster comparison");

  auto mag_mode =
    "mag:" % (clipp::command("mag").set(selected, mode::mag),
              (clipp::command("date") | clipp::command("content")),
              clipp::option("-b", "--binary") % "compare files byte by byte",
              clipp::option("-q", "--quick") % "use heuristics for faster comparison");

  bool version = false, help = false;
  auto cli = (dem_mode | ele_mode | mag_mode,
              clipp::option("-v", "--version").set(version) % "show version",
              clipp::option("-h", "--help").set(help) % "show help");

  auto call_help = [&] {
    auto fmt = clipp::doc_formatting{}.doc_column(31);
    std::cout << make_man_page(cli, argv[0], fmt) << '\n';
  };

  if (clipp::parse(argc, argv, cli))
  {
    switch (selected)
    {
      case mode::dem:
        std::cout << "dem: not implemented" << '\n';
        break;
      case mode::ele:
        std::cout << "ele: not implemented" << '\n';
        break;
      case mode::mag:
        std::cout << "mag: not implemented" << '\n';
        break;
      case mode::help:
        std::cout << "mag: not implemented" << '\n';
        break;
    }
    if (version) std::cout << "version 0.1.0\n\n";
    if (help) call_help();
  }
  else
  {
    call_help();
  }

  // argparse::ArgumentParser program("emgsim");
  // program.add_argument("{dem,ele,mag}").help("sub-command. see each detail '** -h'");
  // program.add_argument("").remaining();
  // program.parse(argc, argv);

  // if (auto sub_command = program.get("{dem,ele,mag}"); sub_command == "dem")
  // {
  //   argparse::ArgumentParser command("emgsim dem");
  //   command.add_argument("-i","--input").help("input file");
  //   auto [sub_argc, sub_argv] = program.get_sub_arg();
  //   std::cout << sub_argc << '\n';
  //   std::cout << sub_argv[0] << '\n';
  //   // for (size_t i = 0; i < sub_argc; ++i)
  //   // {
  //   //   std::cout << sub_argv[i] << '\n';
  //   // }

  //   // command.parse(sub_argc, sub_argv);
  // }
  // else if (sub_command == "ele")
  // {
  //   std::cout << "ele" << '\n';
  // }
  // else if (sub_command == "mag")
  // {
  //   std::cout << "mag" << '\n';
  // }
  // else
  // {
  //   std::cout << ": wrong Positional argument. expected: {dem,ele,mag}" << '\n';
  //   std::cout << program;
  // }

  // dem::DEM daf;
  // daf.SetPathInput(argc, argv);  // load command line input
  // daf.LoadParamater();  // load ini input file
  // daf.SetPathOutput();  // create output directories
  // daf.LoadExtData();  // load magnetic and electrostatic field, particle,
  //                     // object, imaginary object
  // daf.AssociateObject();  // associate same instance information
  // daf.CreateCell();  // create virtual cell for efficient calculation of
  //                    // particle collision, dipole interactions, coulomb force
  //                    // between particles
  // daf.RunDEM();  // dem start

  return 0;
}
