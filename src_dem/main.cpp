#include "DEM.h"

int main(int argc, char** argv) {

	std::ios::sync_with_stdio(false);

	dem::DEM daf;

	daf.SetPathInput(argc, argv);	//load command line input

	daf.LoadParamater();			//load ini input file

	daf.SetPathOutput();			//create output directories

	daf.LoadExtData();				//load magnetic and electrostatic field, particle, object, imaginary object

	daf.AssociateObject();			//associate same instance information

	daf.CreateCell();				//create virtual cell for efficient calculation of particle collision, dipole interactions, coulomb force between particles

	daf.RunDEM();					//dem start

	return 0;
}
