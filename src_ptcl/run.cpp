#include "run.h"

Run::Run() :PI(3.14159265358979323846264338328)
{
	min_Pos[0] = 0;
	max_Pos[0] = 0;
	min_Pos[1] = 0;
	max_Pos[1] = 0;
	min_Pos[2] = 0;
	max_Pos[2] = 0;
	ptclNum = 0;
	load_ptcl_number_diameter = false;
}

Run::~Run()
{

}

void Run::Load(char *filename_1, char *filename_2){

	string temp[50], split;
	int r = 0, p = 0;

	cout << "?¿?v?¿?Z?¿½?¿½?¿½?¿½?¿½?¿½?????¿½?¿½?¿½Å??¿½?¿½???¿?..." << endl;
	ifstream ifs(filename_1);

	//?¿?G?¿½?¿½?¿?[?¿½?¿½?¿½?¿½
	if (!ifs){
		cout << "Error.Can't open input.ini data file." << endl;
		exit(1);
	}

	//?¿½Ç????¿½?¿½???¿½?¿½?¿½?¿½??A1?¿?s?¿½Ú????¿½?¿½?¿½?¿½?¿½?¿½?¿?
	getline(ifs, split);

	//1?¿?s?¿½?¿½?¿½Ì?f?¿?[?¿?^?¿½?¿½?¿½Â??????¿½?¿½??Atemp?¿½É?????¿½Û??¿½?¿½?¿½?¿½???¿½?¿½?¿?
	while (!ifs.eof() && getline(ifs, split)) {
		p = (int)split.find("=");
		temp[r] = split.substr(p + 1);
		r++;
	}

	ifs.close();

	load_ptcl_number_diameter = d.bool_string(temp[0]);
	ptclNum = d.int_string(temp[1]);
	min_Pos[0] = d.double_string(temp[2]);
	max_Pos[0] = d.double_string(temp[3]);
	min_Pos[1] = d.double_string(temp[4]);
	max_Pos[1] = d.double_string(temp[5]);
	min_Pos[2] = d.double_string(temp[6]);
	max_Pos[2] = d.double_string(temp[7]);
	min_Velo[0] = d.double_string(temp[8]);
	max_Velo[0] = d.double_string(temp[9]);
	min_Velo[1] = d.double_string(temp[10]);
	max_Velo[1] = d.double_string(temp[11]);
	min_Velo[2] = d.double_string(temp[12]);
	max_Velo[2] = d.double_string(temp[13]);
	min_omega[0] = d.double_string(temp[14]);
	max_omega[0] = d.double_string(temp[15]);
	min_omega[1] = d.double_string(temp[16]);
	max_omega[1] = d.double_string(temp[17]);
	min_omega[2] = d.double_string(temp[18]);
	max_omega[2] = d.double_string(temp[19]);
	ave_sp_gravity = d.double_string(temp[20]);
	standev_sp_gravity = d.double_string(temp[21]);
	diameter_average = d.double_string(temp[22]);
	diameter_stadard_deviation = d.double_string(temp[23]);
	ave_Charge = d.double_string(temp[24]);
	standev_Charge = d.double_string(temp[25]);
	ave_ele_conductivity = d.double_string(temp[26]);
	standev_ele_conductivity = d.double_string(temp[27]);
	ave_permittivity = d.double_string(temp[28]);
	standev_permittivity = d.double_string(temp[29]);
	ave_permeability = d.double_string(temp[30]);
	standev_permeability = d.double_string(temp[31]);
	ave_adhesion = d.double_string(temp[32]);
	standev_adhesion = d.double_string(temp[33]);
	ave_rolling_fri = d.double_string(temp[34]);
	standev_rolling_fri = d.double_string(temp[35]);

	if (load_ptcl_number_diameter){

		//?¿½?¿½?¿?q?¿?f?¿?[?¿?^?¿½Ì?I?¿?[?¿?v?¿½?¿½
		cout << "?¿½?¿½?¿?q?¿?f?¿?[?¿?^?¿½?¿½?????¿½?¿½?¿½Å??¿½?¿½???¿?..." << endl;
		ifstream fin(filename_2);

		//?¿?G?¿½?¿½?¿?[?¿?`?¿?F?¿?b?¿?N
		if (!fin){
			cout << "Error.Can't open 1ptcl data file." << endl;
			exit(1);
		}

		string temp1, temp2;
		int p = 0;

		//1?¿?s?¿½Ú????¿½?¿½l?¿½?¿½?????¿½?¿½??A?¿½?¿½?¿?q?¿½?¿½?¿½?¿½???¿?
		getline(fin, temp1);
		p = (int)(temp1.find(","));
		temp2 = temp1.substr(0, p);
		ptclNum = d.int_string(temp2);

		Pos_X.assign(ptclNum, 0.0);
		Pos_Y.assign(ptclNum, 0.0);
		Pos_Z.assign(ptclNum, 0.0);
		Velo_X.assign(ptclNum, 0.0);
		Velo_Y.assign(ptclNum, 0.0);
		Velo_Z.assign(ptclNum, 0.0);
		Omega_X.assign(ptclNum, 0.0);
		Omega_Y.assign(ptclNum, 0.0);
		Omega_Z.assign(ptclNum, 0.0);
		Diameter.assign(ptclNum, 0.0);
		Diameter.assign(ptclNum, 0.0);
		SP_Gravity.assign(ptclNum, 0.0);
		Charge.assign(ptclNum, 0.0);
		Conductivity.assign(ptclNum, 0.0);
		Permittivity.assign(ptclNum, 0.0);
		Permeability.assign(ptclNum, 0.0);
		Adhesion.assign(ptclNum, 0.0);
		Rolling_Fri.assign(ptclNum, 0.0);

		//2?¿?s?¿½Ú??¿½Ç????¿½?¿½?¿½Å?A?¿½?¿½?¿½?¿½?¿½?¿½?¿½?¿½
		getline(fin, temp1);

		//?¿½È?~?¿?A?¿?e?¿½?¿½?¿?q?¿½Ì?p?¿½?¿½?¿½?¿½?¿?[?¿?^?¿½?¿½???¿?
		for (int i = 0; i < ptclNum; i++){

			//?¿½?¿½s?¿½Ü??¿½Ü??¿½Ç????¿½?¿½?¿?
			getline(fin, temp1);

			//split?¿½?¿½?¿½?¿½?¿?s?¿?@getline?¿½Å??????¿½?¿½?????¿½?¿½?¿½?¿½?¿½?¿½","?¿½Å??¿½Ø??¿?Astr_List?¿½É????¿?
			list<string> str_List = d.split(temp1, ",");

			//?¿?C?¿?e?¿½?¿½?¿?[?¿?^?¿½?¿½?¿½???  ","?¿½Å??¿½Ø??¿½??½?¿½f?¿?[?¿?^?¿½Ì??¿½?¿½??A1?¿½Ô????¿½?¿½???¿½?¿½????C?¿?e?¿½?¿½?¿?[?¿?^?¿½?¿½?¿½?¿½?¿½??¹?¿½?¿?
			list<string>::iterator iter = str_List.begin();

			//?¿?C?¿?e?¿½?¿½?¿?[?¿?^?¿½?¿½i?¿½ß???A?¿½?¿½?¿½Ì?f?¿?[?¿?^?¿½?¿½?¿½?¿½?¿½?¿½?¿½?¿½@?¿½È?~?¿½Í?J?¿½?¿½???¿?
			++iter;
			++iter;
			++iter;
			++iter;
			++iter;
			++iter;
			++iter;
			++iter;
			++iter;
			++iter;

			temp2 = *iter;
			Diameter[i] = d.double_string(temp2);
			++iter;
		}

		fin.close();
		sort(Diameter.begin(),Diameter.end());     // ?¿?\?¿?[?¿?g?¿½?¿½?¿½?¿½

	}
	else{

		Pos_X.assign(ptclNum, 0.0);
		Pos_Y.assign(ptclNum, 0.0);
		Pos_Z.assign(ptclNum, 0.0);
		Velo_X.assign(ptclNum, 0.0);
		Velo_Y.assign(ptclNum, 0.0);
		Velo_Z.assign(ptclNum, 0.0);
		Omega_X.assign(ptclNum, 0.0);
		Omega_Y.assign(ptclNum, 0.0);
		Omega_Z.assign(ptclNum, 0.0);
		Diameter.assign(ptclNum, 0.0);
		Diameter.assign(ptclNum, 0.0);
		SP_Gravity.assign(ptclNum, 0.0);
		Charge.assign(ptclNum, 0.0);
		Conductivity.assign(ptclNum, 0.0);
		Permittivity.assign(ptclNum, 0.0);
		Permeability.assign(ptclNum, 0.0);
		Adhesion.assign(ptclNum, 0.0);
		Rolling_Fri.assign(ptclNum, 0.0);
	}

	cout << " ?¿½?¿½?¿½?¿½?¿½Ý?u?¿½?¿½?¿?q?¿½?¿½ : " << ptclNum << endl;
}

void Run::Particle_Initial_Placement(){

	if (!load_ptcl_number_diameter){
		SizeCal();		//?¿½Å??¿½?¿½???¿½?¿½a?¿?v?¿?Z
	}

	Sp_gravity();
	ChargeCal();
	EleconductCal();
	PermitCal();
	PermeaCal();
	AdhesionCal();
	RollingFriction();

	VeloCal();
	OmegaCal();
	PlaceCal();
}

void Run::Sp_gravity(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(ave_sp_gravity, standev_sp_gravity);
	double temp;

	for (int i = 0; i < ptclNum; i++){

		do{
			temp = rnd(mt);

		} while (temp < 0);

		SP_Gravity[i] = temp;
	}
}

void Run::SizeCal(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(diameter_average, diameter_stadard_deviation);
	double temp;

	for (int i = 0; i < ptclNum; i++){

		do{
			temp = rnd(mt);

		} while (temp < 0);

		Diameter[i] = temp;
	}

	sort(Diameter.begin(), Diameter.end());     // ?¿?\?¿?[?¿?g?¿½?¿½?¿½?¿½
}

void Run::ChargeCal(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(ave_Charge, standev_Charge);

	for (int i = 0; i < ptclNum; i++){
		Charge[i] = rnd(mt);
	}
}

void Run::EleconductCal(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(ave_ele_conductivity, standev_ele_conductivity);
	double temp;

	for (int i = 0; i < ptclNum; i++){

		do{
			temp = rnd(mt);

		} while (temp < 0);

		Conductivity[i] = temp;
	}
}
void Run::PermitCal(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(ave_permittivity, standev_permittivity);
	double temp;

	for (int i = 0; i < ptclNum; i++){

		do{
			temp = rnd(mt);

		} while (temp < 0);

		Permittivity[i] = temp;
	}
}

void Run::PermeaCal(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(ave_permeability, standev_permeability);
	double temp;

	for (int i = 0; i < ptclNum; i++){

		do{
			temp = rnd(mt);

		} while (temp < 0);

		Permeability[i] = temp;
	}
}

void Run::AdhesionCal(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(ave_adhesion, standev_adhesion);
	double temp;

	for (int i = 0; i < ptclNum; i++){

		do{
			temp = rnd(mt);

		} while (temp < 0);

		Adhesion[i] = temp;
	}

}

void Run::RollingFriction(){

	random_device rnddev;
	mt19937 mt(rnddev());
	normal_distribution<double> rnd(ave_rolling_fri, standev_rolling_fri);
	double temp;

	for (int i = 0; i < ptclNum; i++){

		do{
			temp = rnd(mt);

		} while (temp < 0);

		Rolling_Fri[i] = temp;
	}
}

void Run::PlaceCal(){

	//?¿½?¿½?¿?q?¿½?¿½?¿½?¿½?¿½Ý?u?¿½????¿½?¿½?¿½?¿½
	random_device rnddev;
	mt19937 mt(rnddev());
	uniform_real_distribution<double> rnd_length(min_Pos[0], max_Pos[0]), rnd_depth(min_Pos[1], max_Pos[1]), rnd_height(min_Pos[2], max_Pos[2]);
	double pos_x, pos_y, pos_z;
	const double radius_sphere = (0.025 - 0.0008) * (0.025 - 0.0008);
	const double center_sphere_x = 0.03;
	const double center_sphere_y = 0.03;
	const double center_sphere_z = 0.03;

	for (int i = 0; i < ptclNum; i++){

		temp :
		pos_x = rnd_length(mt);
		pos_y = rnd_depth(mt);
		pos_z = rnd_height(mt);

		if ((pos_x - center_sphere_x) * (pos_x - center_sphere_x) + (pos_y - center_sphere_y) * (pos_y - center_sphere_y) + (pos_z - center_sphere_z) * (pos_z - center_sphere_z) > radius_sphere) goto temp;

		for (int j = i - 1; j >= 0; j--){

			//?¿½?¿½?¿½Î???u
			const double dx = Pos_X[j] - pos_x;
			const double dy = Pos_Y[j] - pos_y;
			const double dz = Pos_Z[j] - pos_z;

			if ((dx*dx + dy*dy + dz*dz) < ((Diameter[i] + Diameter[j])*0.5)*((Diameter[i] + Diameter[j])*0.5)) goto temp;
		}

		Pos_X[i] = pos_x;
		Pos_Y[i] = pos_y;
		Pos_Z[i] = pos_z;

		cout << "?¿?c?¿½??±?¿½q?¿?z?¿?u?¿½?¿½ = " << ptclNum - i << '\r';
	}
}

void Run::VeloCal(){

	//?¿½?¿½?¿?q?¿½?¿½?¿½?¿½?¿½?¿½?¿?x?¿½?¿½?¿½?¿½
	random_device rnddev;
	mt19937 mt(rnddev());
	uniform_real_distribution<double> rnd_length(min_Velo[0], max_Velo[0]), rnd_depth(min_Velo[1], max_Velo[1]), rnd_height(min_Velo[2], max_Velo[2]);

	for (int i = 0; i < ptclNum; i++){

		Velo_X[i] = rnd_length(mt);
		Velo_Y[i] = rnd_depth(mt);
		Velo_Z[i] = rnd_height(mt);
	}
}

void Run::OmegaCal(){

	random_device rnddev;
	mt19937 mt(rnddev());
	uniform_real_distribution<double> rnd_length(min_omega[0], max_omega[0]), rnd_depth(min_omega[1], max_omega[1]), rnd_height(min_omega[2], max_omega[2]);

	for (int i = 0; i < ptclNum; i++){

		Omega_X[i] = rnd_length(mt);
		Omega_Y[i] = rnd_depth(mt);
		Omega_Z[i] = rnd_height(mt);
	}
}

void Run::Output_ptcl(char * filename){

	ofstream ofs;

	ofs.open(filename, ios_base::out, ios_base::trunc);
	ofs << ptclNum << "," << "number" << endl;
	ofs << "x (m),y (m),z (m),u (m/s),v (m/s),w (m/s),omega1 (rad/s),omega2 (rad/s),omega3 (rad/s),sp_gravity (kg/m3),diameter (m),charge/mass (?¿½?¿½C/g),conductivity (?¿½?¿½-1m-1),permittivity,permeability,adhesion,rolling_friction" << endl;

	for (int i = 0; i<ptclNum; i++){

		ofs << Pos_X[i] << "," << Pos_Y[i] << "," << Pos_Z[i] << "," << Velo_X[i] << "," << Velo_Y[i] << "," << Velo_Z[i]
			<< "," << Omega_X[i] << "," << Omega_Y[i] << "," << Omega_Y[i] << "," << SP_Gravity[i] << "," << Diameter[i]
			<< "," << Charge[i] << "," << Conductivity[i] << "," << Permittivity[i] << "," << Permeability[i]
			<< "," << Adhesion[i] << "," << Rolling_Fri[i] << endl;
	}

	ofs.close();
}
