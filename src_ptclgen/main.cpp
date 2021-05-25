#include "Generator.h"

int main(int argc, char** argv) {

	std::ios::sync_with_stdio(false);
	dem::Generator gen;
  gen.read_args(argc, argv);
  gen.set();
  gen.generate();
	return 0;
}
