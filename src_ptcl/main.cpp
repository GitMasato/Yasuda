#include "run.h"

int main(int argc, char **argv)
{
  std::ios::sync_with_stdio(false);

  Run run;
	char input_file_1[64] = ".\\INPUT\\Load.ini";
	char input_file_2[64] = ".\\INPUT\\1ptcl.csv";
	char out_putfile[64] = ".\\OUTPUT\\1ptcl.csv";
	run.Load(input_file_1, input_file_2);
	run.Particle_Initial_Placement();
	run.Output_ptcl(out_putfile);
}
