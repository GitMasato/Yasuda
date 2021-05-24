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

	cout << "計算条件を読み込んでいます..." << endl;
	ifstream ifs(filename_1);

	//エラー処理
	if (!ifs){
		cout << "Error.Can't open input.ini data file." << endl;
		exit(1);
	}

	//読み込むだけで、1行目は無視する
	getline(ifs, split);

	//1行分のデータずつ読み込み、tempに一時保存していく
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

		//粒子データのオープン
		cout << "粒子データを読み込んでいます..." << endl;
		ifstream fin(filename_2);

		//エラーチェック
		if (!fin){
			cout << "Error.Can't open 1ptcl data file." << endl;
			exit(1);
		}

		string temp1, temp2;
		int p = 0;

		//1行目の数値を読み込み、粒子数を設定
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

		//2行目を読み込んで、無視する
		getline(fin, temp1);

		//以降、各粒子のパラメータを設定
		for (int i = 0; i < ptclNum; i++){

			//一行まるまる読み込む
			getline(fin, temp1);

			//splitを実行　getlineで読み込んだ文字列を","で区切り、str_Listに保存
			list<string> str_List = d.split(temp1, ",");

			//イテレータを取得  ","で区切られたデータの中で、1番最初のものにイテレータを合わせる
			list<string>::iterator iter = str_List.begin();

			//イテレータを進めて、次のデータを代入する　以降は繰り返し
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
		sort(Diameter.begin(),Diameter.end());     // ソートする

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
	
	cout << " 初期設置粒子数 : " << ptclNum << endl;
}

void Run::Particle_Initial_Placement(){

	if (!load_ptcl_number_diameter){
		SizeCal();		//最初に粒径計算
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
	
	sort(Diameter.begin(), Diameter.end());     // ソートする
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
		
	//粒子初期設置場所決定
	random_device rnddev;
	mt19937 mt(rnddev());
	uniform_real_distribution<double> rnd_length(min_Pos[0], max_Pos[0]), rnd_depth(min_Pos[1], max_Pos[1]), rnd_height(min_Pos[2], max_Pos[2]);
	double pos_x, pos_y, pos_z;

	for (int i = 0; i < ptclNum; i++){				
				
		temp :
		pos_x = rnd_length(mt);
		pos_y = rnd_depth(mt);
		pos_z = rnd_height(mt);

		for (int j = i - 1; j >= 0; j--){			

			//相対位置
			const double dx = Pos_X[j] - pos_x;
			const double dy = Pos_Y[j] - pos_y;
			const double dz = Pos_Z[j] - pos_z;

			if ((dx*dx + dy*dy + dz*dz) < ((Diameter[i] + Diameter[j])*0.5)*((Diameter[i] + Diameter[j])*0.5)) goto temp;
		}
		
		Pos_X[i] = pos_x;
		Pos_Y[i] = pos_y;
		Pos_Z[i] = pos_z;

		cout << "残り粒子配置数 = " << ptclNum - i << '\r';
	}
}

void Run::VeloCal(){

	//粒子初期速度決定
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

	//粒子初期角速度決定
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
	ofs << "x (m),y (m),z (m),u (m/s),v (m/s),w (m/s),omega1 (rad/s),omega2 (rad/s),omega3 (rad/s),sp_gravity (kg/m3),diameter (m),charge/mass (μC/g),conductivity (Ω-1m-1),permittivity,permeability,adhesion,rolling_friction" << endl;

	for (int i = 0; i<ptclNum; i++){

		ofs << Pos_X[i] << "," << Pos_Y[i] << "," << Pos_Z[i] << "," << Velo_X[i] << "," << Velo_Y[i] << "," << Velo_Z[i]
			<< "," << Omega_X[i] << "," << Omega_Y[i] << "," << Omega_Y[i] << "," << SP_Gravity[i] << "," << Diameter[i]
			<< "," << Charge[i] << "," << Conductivity[i] << "," << Permittivity[i] << "," << Permeability[i] 
			<< "," << Adhesion[i] << "," << Rolling_Fri[i] << endl;
	}
	
	ofs.close();
}