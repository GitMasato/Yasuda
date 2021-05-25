#include "Ele_Calc.h"

Electric_Calc::Electric_Calc() :PI(3.14159265358979323846264338328){		//コンストラクタ．変数の初期化をしている．

	//基本パラメータ
	ElmSize = 200.0e-6;			//メッシュサイズ
	areaWidth = 50000.0e-6;		//Ｘ方向計算領域
	areaDepth = 50000.0e-6;		//Y方向の計算領域
	areaHeight = 50000.0e-6;	//Ｚ方向計算領域	
	Type_boundary_X = 0;		//隅の境界条件の設定　0 = ノイマン / 1 = 周期 から選択
	Type_boundary_Y = 0;		//隅の境界条件の設定　0 = ノイマン / 1 = 周期 から選択

	spElePerm = 8.854187816e-12;//誘電率
	AirRelElePerm = 1.00059;	//空間の比誘電率
	Volt = 5000;				//印加電圧		
	TFC = 1.0e-3;				//収束判定
	SOR = 1.8;					//加速係数	
	OMP_CALC = true;			//並列化するか
	thread = 6;					//並列数
	unsteady_calc = false;		//非定常の計算か
	load_ptcl_Poisson = false;	//粒子データを読み込むか　ポアソンか
	
	//非定常電界計算を行う場合のみ
	if (unsteady_calc){
		total_T = 0.5;					// シミュレーション時間								
		delta_T = 0.5e-3;				// 離散化時間は時定数の 1/10 以下になるように
		// (ParticleRelElePerm * spElePerm)/ conductivity)		-6 : 2.21e-5, -7 : 2.21e-4, -8 : 2.21e-3, -9 : 2.21e-2, -10 : 2.21e-1,	-11 : 2.21,   -12 : 2.21e+1		時定数
		frequency = 1;					//周波数
		phase = 2;						//電位の相数
		rise_time = 2.5e-3;				//電位の立ち上がり時間
		WriteAVS = 50;					//動画の書き込みステップ
		total_step = (int)(total_T / delta_T);
		phase_No = 2;
		phase_change_flag = true;
		span = 0.0;
	}

	//粒子を考慮する・ポアソンの場合のみ
	if (load_ptcl_Poisson){
		ParticleRelElePerm = 2.5;		//粒子の比誘電率
		conductivity = 1.0e-9;			//導電率
	}	
	
	//線電極を入れる場合
	cylinder_Xdir_Num = 6 ;					//電極数
	cylinder_Xdir.resize(cylinder_Xdir_Num);		
	CylinderElePerm = 5.0;				//電極被服の誘電率

	for (int i = 0; i < cylinder_Xdir_Num; i++){

		cylinder_Xdir[i].Pos_X = 7500.0e-6;
		cylinder_Xdir[i].Pos_Y = 12500.0e-6 + 5000.0e-6 * i;
		cylinder_Xdir[i].Pos_Z = 9800.0e-6;
		cylinder_Xdir[i].pos_X = (int)(cylinder_Xdir[i].Pos_X / ElmSize + 0.5);
		cylinder_Xdir[i].pos_Y = (int)(cylinder_Xdir[i].Pos_Y / ElmSize + 0.5);
		cylinder_Xdir[i].pos_Z = (int)(cylinder_Xdir[i].Pos_Z / ElmSize + 0.5);

		cylinder_Xdir[i].Radius = 400.0e-6;					//電極半径
		cylinder_Xdir[i].Length = 35000.0e-6;				//電極長さ
		cylinder_Xdir[i].Dielectric_Radius = 900.0e-6;		//電極被服を含めた半径
		cylinder_Xdir[i].radius = (int)(cylinder_Xdir[i].Radius / ElmSize + 0.5);
		cylinder_Xdir[i].length = (int)(cylinder_Xdir[i].Length / ElmSize + 0.5);
		cylinder_Xdir[i].diele_radius = (int)(cylinder_Xdir[i].Dielectric_Radius / ElmSize + 0.5);
	}

	//線電極を入れる場合
	cylinder_Ydir_Num = 6;					//電極数
	cylinder_Ydir.resize(cylinder_Ydir_Num);

	for (int i = 0; i < cylinder_Ydir_Num; i++){

		cylinder_Ydir[i].Pos_X = 12500.0e-6 + 5000.0e-6 * i;
		cylinder_Ydir[i].Pos_Y = 7500.0e-6;
		cylinder_Ydir[i].Pos_Z = 6800.0e-6;
		cylinder_Ydir[i].pos_X = (int)(cylinder_Ydir[i].Pos_X / ElmSize + 0.5);
		cylinder_Ydir[i].pos_Y = (int)(cylinder_Ydir[i].Pos_Y / ElmSize + 0.5);
		cylinder_Ydir[i].pos_Z = (int)(cylinder_Ydir[i].Pos_Z / ElmSize + 0.5);

		cylinder_Ydir[i].Radius = 400.0e-6;				//電極半径
		cylinder_Ydir[i].Length = 35000.0e-6;			//電極長さ
		cylinder_Ydir[i].Dielectric_Radius = 900.0e-6;  //電極被服を含めた半径
		cylinder_Ydir[i].radius = (int)(cylinder_Ydir[i].Radius / ElmSize + 0.5);
		cylinder_Ydir[i].length = (int)(cylinder_Ydir[i].Length / ElmSize + 0.5);
		cylinder_Ydir[i].diele_radius = (int)(cylinder_Ydir[i].Dielectric_Radius / ElmSize + 0.5);
	}	

	//円筒を入れる場合
	hollow_cylNum = 1;	
	hollow_cylinder.resize(hollow_cylNum);
	hollow_cylinderElePerm = 3.1;

	for (int i = 0; i < hollow_cylNum; i++){

		hollow_cylinder[i].Pos_X = 25000.0e-6;
		hollow_cylinder[i].Pos_Y = 25000.0e-6;
		hollow_cylinder[i].Pos_Z = 0.0e-6;
		hollow_cylinder[i].pos_X = (int)(hollow_cylinder[i].Pos_X / ElmSize + 0.5);
		hollow_cylinder[i].pos_Y = (int)(hollow_cylinder[i].Pos_Y / ElmSize + 0.5);
		hollow_cylinder[i].pos_Z = (int)(hollow_cylinder[i].Pos_Z / ElmSize + 0.5);

		hollow_cylinder[i].Inside_Radius = 17500.0e-6;
		hollow_cylinder[i].Outside_Radius = 25000.0e-6;
		hollow_cylinder[i].Length = 16000.0e-6;				
		hollow_cylinder[i].inside_radius = (int)(hollow_cylinder[i].Inside_Radius / ElmSize + 0.5);
		hollow_cylinder[i].outside_radius = (int)(hollow_cylinder[i].Outside_Radius / ElmSize + 0.5);
		hollow_cylinder[i].length = (int)(hollow_cylinder[i].Length / ElmSize + 0.5);
	}

	//砂場を入れる場合
	wallNum = 1;
	wall.resize(wallNum);
	wallElePerm = 2.5;

	for (int i = 0; i < wallNum; i++){

		wall[i].Pos_X[0] = 0;
		wall[i].Pos_Y[0] = 0;
		wall[i].Pos_Z[0] = 0;

		wall[i].Pos_X[1] = areaWidth;
		wall[i].Pos_Y[1] = 0;
		wall[i].Pos_Z[1] = 0;

		wall[i].Pos_X[2] = 0;
		wall[i].Pos_Y[2] = 0;
		wall[i].Pos_Z[2] = 5900.0e-6;

		wall[i].Pos_X[3] = areaWidth;
		wall[i].Pos_Y[3] = 0;
		wall[i].Pos_Z[3] = 5900.0e-6;		

		for (int n = 0; n < 4; n++){
			wall[i].pos_X[n] = (int)(wall[i].Pos_X[n] / ElmSize + 0.5);
			wall[i].pos_Y[n] = (int)(wall[i].Pos_Y[n] / ElmSize + 0.5);
			wall[i].pos_Z[n] = (int)(wall[i].Pos_Z[n] / ElmSize + 0.5);
		}

		wall[i].Length = areaDepth;
		wall[i].length = (int)(wall[i].Length / ElmSize + 0.5);		
	}

	//以降は変数の初期化・計算パラメータの離散化
	//基本パラメータに関して
	WidthNum = (int)(areaWidth / ElmSize + 0.5);
	DepthNum = (int)(areaDepth / ElmSize + 0.5);
	HeightNum = (int)(areaHeight / ElmSize + 0.5);
	cellNum_2D = WidthNum * HeightNum;
	cellNum_3D = WidthNum * DepthNum * HeightNum;
	if (cellNum_3D == 0)	cellNum_3D = 1;
	convergentCount = 0;

	ElePerm.assign(cellNum_3D, 0.0);	Potential.assign(cellNum_3D, 0.0);
	E_x.assign(cellNum_3D, 0.0);	E_y.assign(cellNum_3D, 0.0);	E_z.assign(cellNum_3D, 0.0);
	bcc.assign(cellNum_3D, 0);		fcc.assign(cellNum_3D, 0);		cbc.assign(cellNum_3D, 0);
	cfc.assign(cellNum_3D, 0);		ccb.assign(cellNum_3D, 0);		ccf.assign(cellNum_3D, 0);
	bbb.assign(cellNum_3D, 0);		bbf.assign(cellNum_3D, 0);		bfb.assign(cellNum_3D, 0);		bff.assign(cellNum_3D, 0);
	fbb.assign(cellNum_3D, 0);		fbf.assign(cellNum_3D, 0);		ffb.assign(cellNum_3D, 0);		fff.assign(cellNum_3D, 0);

	CalcFlag = new bool[cellNum_3D];
	for (int i = 0; i < cellNum_3D; i++)	CalcFlag[i] = true;

	//粒子を考慮する場合・ポアソンの式の場合
	if (load_ptcl_Poisson){
		ptclNum = 0;
		ptcl_num_node.assign(cellNum_3D, -1);
		ptcl_num_ele.assign(cellNum_3D, -1);
		Conductivity.assign(cellNum_3D, 0);
		Charge.assign(cellNum_3D, 0);
	}

	cout << "計算条件を設定完了" << endl;
}

Electric_Calc::~Electric_Calc(){

	if (CalcFlag !=0)	delete[] CalcFlag;
}

void Electric_Calc::Setting_Boundary_Condition_V(){

	//隅の電位に関する境界条件設定
	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				bcc[m] = m - 1;
				fcc[m] = m + 1;
				cbc[m] = m - WidthNum;
				cfc[m] = m + WidthNum;
				ccb[m] = m - (WidthNum * DepthNum);
				ccf[m] = m + (WidthNum * DepthNum);
			}
		}
	}

	if (Type_boundary_X == 0)			Setting_Boundary_Condition_V_X_n();
	else if (Type_boundary_X == 1)		Setting_Boundary_Condition_V_X_l();
	
	if (Type_boundary_Y == 0)			Setting_Boundary_Condition_V_Y_n();
	else if (Type_boundary_Y == 1)		Setting_Boundary_Condition_V_Y_l();
	
	Setting_Boundary_Condition_V_Z_n();
	
	cout << "電位に関する端の境界条件設定" << endl;
}

void Electric_Calc::Setting_Boundary_Condition_V_X_n(){

	for(int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (i == 0)					bcc[m] = m + 1;
				if (i == WidthNum - 1)		fcc[m] = m - 1;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_X_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (i == 0)					bcc[m] = m + WidthNum - 1;
				if (i == WidthNum - 1)		fcc[m] = m - WidthNum + 1;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_Y_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (j == 0)					cbc[m] = m + WidthNum;
				if (j == DepthNum - 1)		cfc[m] = m - WidthNum;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_Y_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (j == 0)					cbc[m] = m + (HeightNum * WidthNum) - WidthNum;
				if (j == DepthNum - 1)		cfc[m] = m - (HeightNum * WidthNum) + WidthNum;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (k == 0)					ccb[m] = m + (WidthNum * DepthNum);
				if (k == HeightNum - 1)		ccf[m] = m - (WidthNum * DepthNum);
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e(){

	//隅の誘電率に関する境界条件設定
	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				bbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
				fbb[m] = m - (WidthNum*DepthNum) - WidthNum;
				bfb[m] = m - (WidthNum*DepthNum) - 1;
				ffb[m] = m - (WidthNum*DepthNum);
				bbf[m] = m - WidthNum - 1;
				fbf[m] = m - WidthNum;
				bff[m] = m - 1;
				fff[m] = m;
			}
		}
	}
	Setting_Boundary_Condition_e_X_n_corner();          //以下端と角の誘電率に関しては計算に大きな影響を与えないため、ノイマンとして扱う(電位に関しては考慮すべき)
	Setting_Boundary_Condition_e_Y_n_corner();          //以下端と角の誘電率に関しては計算に大きな影響を与えないため、ノイマンとして扱う
	Setting_Boundary_Condition_e_Z_n_corner();          //以下端と角の誘電率に関しては計算に大きな影響を与えないため、ノイマンとして扱う
	Setting_Boundary_Condition_e_X_n_Y_n_corner();      //以下端と角の誘電率に関しては計算に大きな影響を与えないため、ノイマンとして扱う
	Setting_Boundary_Condition_e_X_n_Z_n_corner();      //以下端と角の誘電率に関しては計算に大きな影響を与えないため、ノイマンとして扱う
	Setting_Boundary_Condition_e_Y_n_Z_n_corner();      //以下端と角の誘電率に関しては計算に大きな影響を与えないため、ノイマンとして扱う
	Setting_Boundary_Condition_e_X_n_Y_n_Z_n_corner(); 

	cout << "誘電率に関する端の境界条件設定" << endl;
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (i == 0){
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bfb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - WidthNum;
					bff[m] = m;
				}
				if (i == WidthNum - 1){
					fbb[m] = m - (WidthNum * DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum * DepthNum) - 1;
					fbf[m] = m - WidthNum - 1;
					fff[m] = m - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_Y_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (j == 0){
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - 1;
					fbf[m] = m;
				}
				if (j == DepthNum - 1){
					bfb[m] = m - (WidthNum * DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum * DepthNum) - WidthNum;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (k == 0){
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum;
					bfb[m] = m - 1;
					ffb[m] = m;
				}
				if (k == HeightNum - 1){
					bbf[m] = m - (WidthNum * DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum * DepthNum) - WidthNum;
					bff[m] = m - (WidthNum * DepthNum) - 1;
					fff[m] = m - (WidthNum * DepthNum);
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_Y_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((i == 0) && (j == 0)){
					bbb[m] = m - (WidthNum*DepthNum);
					fbb[m] = m - (WidthNum*DepthNum);
					bfb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m;
					fbf[m] = m;
					bff[m] = m;
				}
				if ((i == 0) && (j == DepthNum - 1)){
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bbf[m] = m - WidthNum;
					bff[m] = m - WidthNum;
					fff[m] = m - WidthNum;
				}
				if ((i == WidthNum - 1) && (j == 0)){
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum) - 1;
					ffb[m] = m - (WidthNum*DepthNum) - 1;
					bbf[m] = m - 1;
					fbf[m] = m - 1;
					fff[m] = m - 1;
				}
				if ((i == WidthNum - 1) && (j == DepthNum - 1)){
					fbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - WidthNum - 1;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((i == 0) && (k == 0)){
					bbb[m] = m - WidthNum;
					fbb[m] = m - WidthNum;
					bfb[m] = m;
					ffb[m] = m;
					bbf[m] = m - WidthNum;					
					bff[m] = m;
				}
				if ((i == 0) && (k == HeightNum - 1)){
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;					
					bfb[m] = m - (WidthNum*DepthNum);					
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					bff[m] = m - (WidthNum*DepthNum);
					fff[m] = m - (WidthNum*DepthNum);
				}
				if ((i == WidthNum - 1) && (k == 0)){
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum - 1;
					bfb[m] = m - 1;
					ffb[m] = m - 1;					
					fbf[m] = m - WidthNum - 1;					
					fff[m] = m - 1;
				}
				if ((i == WidthNum - 1) && (k == HeightNum - 1)){					
					fbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;					
					ffb[m] = m - (WidthNum*DepthNum) - 1;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bff[m] = m - (WidthNum*DepthNum) - 1;
					fff[m] = m - (WidthNum*DepthNum) - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_Y_n_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((j == 0) && (k == 0)){
					bbb[m] = m - 1;
					fbb[m] = m;
					bfb[m] = m - 1;
					ffb[m] = m;
					bbf[m] = m - 1;
					fbf[m] = m;
				}
				if ((j == 0) && (k == HeightNum - 1)){
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - (WidthNum*DepthNum) - 1;
					fbf[m] = m - (WidthNum*DepthNum);
					bff[m] = m - (WidthNum*DepthNum) - 1;
					fff[m] = m - (WidthNum*DepthNum);
				}
				if ((j == DepthNum - 1) && (k == 0)){
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum;
					bfb[m] = m - WidthNum - 1;
					ffb[m] = m - WidthNum;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum;
				}
				if ((j == DepthNum - 1) && (k == HeightNum - 1)){
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					bff[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fff[m] = m - (WidthNum*DepthNum) - WidthNum;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_Y_n_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((i == 0) && (j == 0) && (k == 0)){//fff
					bbb[m] = m;
					fbb[m] = m;
					bfb[m] = m;
					ffb[m] = m;
					bbf[m] = m;
					fbf[m] = m;
					bff[m] = m;
				}
				if ((i == 0) && (j == 0) && (k == HeightNum - 1)){//ffb
					bbb[m] = m - (WidthNum*DepthNum);
					fbb[m] = m - (WidthNum*DepthNum);
					bfb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - (WidthNum*DepthNum);
					fbf[m] = m - (WidthNum*DepthNum);
					bff[m] = m - (WidthNum*DepthNum);
					fff[m] = m - (WidthNum*DepthNum);
				}
				if ((i == 0) && (j == DepthNum -1) && (k == 0)){//fbf
					bbb[m] = m - WidthNum;
					fbb[m] = m - WidthNum;
					bfb[m] = m - WidthNum;
					ffb[m] = m - WidthNum;
					bbf[m] = m - WidthNum;
					bff[m] = m - WidthNum;
					fff[m] = m - WidthNum;
				}
				if ((i == 0) && (j == DepthNum - 1) && (k == HeightNum -1)){//fbb
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					bff[m] = m - (WidthNum*DepthNum) - WidthNum;
					fff[m] = m - (WidthNum*DepthNum) - WidthNum;
				}
				if ((i == WidthNum - 1) && (j == 0) && (k == 0)){//bff
					bbb[m] = m - 1;
					fbb[m] = m - 1;
					bfb[m] = m - 1;
					ffb[m] = m - 1;
					bbf[m] = m - 1;
					fbf[m] = m - 1;
					fff[m] = m - 1;
				}
				if ((i == WidthNum - 1) && (j == 0) && (k == HeightNum - 1)){//bfb
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum) - 1;
					ffb[m] = m - (WidthNum*DepthNum) - 1;
					bbf[m] = m - (WidthNum*DepthNum) - 1;
					fbf[m] = m - (WidthNum*DepthNum) - 1;
					bff[m] = m - (WidthNum*DepthNum) - 1;
					fff[m] = m - (WidthNum*DepthNum) - 1;
				}
				if ((i == WidthNum - 1) && (j == DepthNum - 1) && (k == 0)){//bbf
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum - 1;
					bfb[m] = m - WidthNum - 1;
					ffb[m] = m - WidthNum - 1;
					fbf[m] = m - WidthNum - 1;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum - 1;
				}
				if ((i == WidthNum - 1) && (j == DepthNum - 1) && (k == HeightNum - 1)){//bbb
					fbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bff[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fff[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_1(){

	cout << "電極等の境界条件をセッテイング" << endl;

	//境界条件電位セット
	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				ElePerm[m] = spElePerm * AirRelElePerm;			

				//電極情報を設定
				for (int e = 0; e < cylinder_Xdir_Num; e++){

					double dy = abs(j - cylinder_Xdir[e].pos_Y);
					double dz = abs(k - cylinder_Xdir[e].pos_Z);
					double distance = sqrt((dy * dy) + (dz * dz));

					//電位を設定
					if ((distance <= cylinder_Xdir[e].radius) && (cylinder_Xdir[e].pos_X <= i) && (i <= (cylinder_Xdir[e].pos_X + cylinder_Xdir[e].length))){
						cylinder_Xdir[e].cyl_mesh.push_back(m);
						CalcFlag[m] = false;
						Potential[m] = - Volt;
					}

					//誘電率設定
					if ((distance <= cylinder_Xdir[e].diele_radius) && (cylinder_Xdir[e].pos_X <= i) && (i <= (cylinder_Xdir[e].pos_X + cylinder_Xdir[e].length))){
						ElePerm[m] = spElePerm * CylinderElePerm;
					}
				}

				//電極情報を設定
				for (int e = 0; e < cylinder_Ydir_Num; e++){

					double dx = abs(i - cylinder_Ydir[e].pos_X);
					double dz = abs(k - cylinder_Ydir[e].pos_Z);
					double distance = sqrt((dx * dx) + (dz * dz));

					//電位を設定
					if ((distance <= cylinder_Ydir[e].radius) && (cylinder_Ydir[e].pos_Y <= j) && (j <= (cylinder_Ydir[e].pos_Y + cylinder_Ydir[e].length))){
						cylinder_Ydir[e].cyl_mesh.push_back(m);
						CalcFlag[m] = false;
						Potential[m] = Volt;
					}

					//誘電率設定
					if ((distance <= cylinder_Ydir[e].diele_radius) && (cylinder_Ydir[e].pos_Y <= j) && (j <= (cylinder_Ydir[e].pos_Y + cylinder_Ydir[e].length))){
						ElePerm[m] = spElePerm * CylinderElePerm;
					}
				}

				//中空円筒情報を設定
				for (int e = 0; e < hollow_cylNum; e++){

					double dx = abs(i - hollow_cylinder[e].pos_X);
					double dy = abs(j - hollow_cylinder[e].pos_Y);
					double distance = sqrt((dx * dx) + (dy * dy));

					//誘電率設定
					if ((hollow_cylinder[e].inside_radius <= distance)
						&& (distance <= hollow_cylinder[e].outside_radius) 
						&& (hollow_cylinder[e].pos_Z <= k) 
						&& (k <= (hollow_cylinder[e].pos_Z + hollow_cylinder[e].length))){

						hollow_cylinder[e].hollow_cyl_mesh.push_back(m);
						ElePerm[m] = spElePerm * hollow_cylinderElePerm;
					}
				}

				//砂場情報を設定
				for (int e = 0; e < wallNum; e++){

					//誘電率設定
					if ((i >= wall[e].pos_X[0]) && (i <= wall[e].pos_X[1])
						&& (j >= wall[e].pos_Y[0]) && (j <= wall[e].pos_Y[0] + wall[e].length)
						&& (k >= wall[e].pos_Z[0]) && (k <= wall[e].pos_Z[2])){

						wall[e].wall_mesh.push_back(m);
						ElePerm[m] = spElePerm * wallElePerm;
					}
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_2(int step, char *filename){

	//境界条件電位セット
	//電位立ち上がり時
	if (((step * delta_T) >= span) && ((step * delta_T) <= span + rise_time)){

		//電界を切り替えるか
		if (phase_change_flag){

			if (phase_No == 1)			phase_No = 2;
			else if (phase_No == 2)		phase_No = 1;			
			phase_change_flag = false;
		}

		if (phase_No == 1){			
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = -Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}		
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}
		}
		else if (phase_No == 2){			
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = Potential[m] = Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = -Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}
		}

		//電位立ち上がり時間を超えた瞬間
	}
	else if ((step * delta_T) > span + rise_time){

		phase_change_flag = true;
		span += (1.0 / frequency) / phase;		

		if (phase_No == 1){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = -Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = Volt;
				}
			}
		}
		else if (phase_No == 2){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = -Volt;
				}
			}
		}

		//電位立ち上がり時間を超えた後
	}
	else{
		if (phase_No == 1){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = -Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = Volt;
				}
			}
		}
		else if (phase_No == 2){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = -Volt;
				}
			}
		}
	}

	//電位データを出力	
	ofstream fout;
	if (step == 0){
		fout.open(filename, ios_base::out, ios_base::trunc);
	}
	else{
		fout.open(filename, ios_base::app);
	}

	for (int e = 0; e < cylinder_Xdir_Num; e++){
		fout << Potential[cylinder_Xdir[e].cyl_mesh[0]] << ",";
	}

	for (int e = 0; e < cylinder_Ydir_Num; e++){
		fout << Potential[cylinder_Ydir[e].cyl_mesh[0]] << ",";
	}
	fout << endl;
	fout.close();
}

void Electric_Calc::Ptcl_Data(char *filename){

	//粒子データのオープン
	cout << "粒子データを読み込んでいます..." << endl;
	ifstream fin(filename);

	//エラーチェック
	if (!fin){
		cout << "Error.Can't open 1ptcl data file." << endl;
		exit(1);
	}

	string temp1, temp2;
	int p = 0;

	//粒子データの読みこみ
	//1行目の数値を読み込み、粒子数を設定
	getline(fin, temp1);
	p = int(temp1.find(","));
	temp2 = temp1.substr(0, p);
	ptclNum = int_string(temp2);

	P_C.resize(ptclNum);

	for (int i = 0; i < ptclNum; i++){
		P_C[i].Ptcl_Pos_X = 0.0;
		P_C[i].Ptcl_Pos_Y = 0.0;
		P_C[i].Ptcl_Pos_Z = 0.0;
		P_C[i].Ptcl_Radius = 0.0;
		P_C[i].Ptcl_Charge = 0.0;
		P_C[i].sp_G = 0.0;
	}	
	
	//2行目を読み込んで、無視する
	getline(fin, temp1);

	//各粒子のパラメータを設定
	for (int i = 0; i < ptclNum; i++){

		//一行まるまる読み込む
		getline(fin, temp1);

		//splitを実行　getlineで読み込んだ文字列を","で区切り、str_Listに保存
		list<string> str_List = split(temp1, ",");

		//イテレータを取得  ","で区切られたデータの中で、1番最初のものにイテレータを合わせる
		list<string>::iterator iter = str_List.begin();

		//イテレータが指示するデータを代入する
		temp2 = *iter;
		P_C[i].Ptcl_Pos_X = double_string(temp2);
		++iter;

		//イテレータを進めて、次のデータを代入する　以降は繰り返し
		temp2 = *iter;
		P_C[i].Ptcl_Pos_Y = double_string(temp2);
		++iter;

		temp2 = *iter;
		P_C[i].Ptcl_Pos_Z = double_string(temp2);
		++iter;

		++iter;
		++iter;
		++iter;
		++iter;
		++iter;
		++iter;

		temp2 = *iter;
		P_C[i].sp_G = double_string(temp2);
		++iter;

		temp2 = *iter;
		P_C[i].Ptcl_Radius = 0.5 * double_string(temp2);
		++iter;

		temp2 = *iter;
		P_C[i].Ptcl_Charge = 4.0 / 3.0 * PI * pow(P_C[i].Ptcl_Radius, 3) * P_C[i].sp_G * double_string(temp2) * 1.0E-3;
	}

	fin.close();
}

void Electric_Calc::Particle_Node_Mapping(){

	cout << "節点上に粒子を設置します" << endl;

	int P_x, P_y, P_z, P_r, m;
	double rel_x, rel_y, rel_z, rel;

	//粒子中心位置を特定
	for (int p = 0; p < ptclNum; p++){

		//粒子位置を整数座標に
		P_x = (int)(P_C[p].Ptcl_Pos_X / ElmSize + 0.5);
		P_y = (int)(P_C[p].Ptcl_Pos_Y / ElmSize + 0.5);
		P_z = (int)(P_C[p].Ptcl_Pos_Z / ElmSize + 0.5);
		P_r = (int)(P_C[p].Ptcl_Radius / ElmSize + 0.5);

		for (int k = 0; k <= HeightNum; k++){
			for (int j = P_y - P_r - 1; j <= P_y + P_r + 1; j++){
				for (int i = P_x - P_r - 1; i <= P_x + P_r + 1; i++){

					//粒子と対象セルの相対距離を算出
					rel_x = (P_x - i) * ElmSize;
					rel_y = (P_y - j) * ElmSize;
					rel_z = (P_z - k) * ElmSize;
					rel = rel_x*rel_x + rel_y*rel_y + rel_z*rel_z;

					//重なっていれば粒子情報をセット
					if (rel == 0){
						m = k * WidthNum * DepthNum + j * WidthNum + i;
						ptcl_num_node[m] = p;
						P_C[p].node.push_back(m);
					}
					else{
						if (rel <= pow(P_C[p].Ptcl_Radius, 2)){

							//節点の番号を求める
							//X方向が-の節点
							if (i<0){
								//Y方向が-の節点
								if (j<0){
									m = k * WidthNum * DepthNum + j * WidthNum + i + WidthNum * DepthNum + WidthNum;
									//Y方向が最大値を超えた節点
								}
								else if (j >= DepthNum){
									m = k * WidthNum * DepthNum + j * WidthNum + i - WidthNum * DepthNum + WidthNum;
								}
								else{
									m = k * WidthNum * DepthNum + j * WidthNum + i + WidthNum;
								}
								//X方向が最大値を超えた節点
							}
							else if (i >= WidthNum){
								//Y方向が-の節点
								if (j<0){
									m = k * WidthNum * DepthNum + j * WidthNum + i + WidthNum * DepthNum - WidthNum;
									//Y方向が最大値を超えた節点
								}
								else if (j >= DepthNum){
									m = k * WidthNum * DepthNum + j * WidthNum + i - WidthNum * DepthNum - WidthNum;
								}
								else{
									m = k * WidthNum * DepthNum + j * WidthNum + i - WidthNum;
								}
							}
							else{
								//Y方向が-の節点
								if (j<0){
									m = k * WidthNum * DepthNum + j * WidthNum + WidthNum * DepthNum + i;
									//Y方向が最大値を超えた節点
								}
								else if (j >= DepthNum){
									m = k * WidthNum * DepthNum + j * WidthNum - WidthNum * DepthNum + i;
								}
								else{
									m = k * WidthNum * DepthNum + j * WidthNum + i;
								}
							}

							ptcl_num_node[m] = p;
							P_C[p].node.push_back(m);
						}
					}
				}
			}
		}
	}

	/*
	
	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open("check.fld", ios_base::out | ios_base::trunc);
	fout_charge2.open("check.dat", ios_base::out | ios_base::trunc);

	fout_charge1
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << "check.dat" << ' ' << "filetype=ascii\n";
	
	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_charge2 << ptcl_num_node[m] << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();

	*/

}

void Electric_Calc::Particle_Ele_Mapping(){

	cout << "要素上に粒子を設置します" << endl;

	int p_node_0, p_node_1, p_node_2, p_node_3, p_node_4, p_node_5, p_node_6, p_node_7;

	//粒子中心位置を特定
	for (int p = 0; p < ptclNum; p++){
		int node = (int)P_C[p].node.size();

		for (int i = 0; i < node; i++){
			int m = P_C[p].node[i];

			p_node_0 = ptcl_num_node[m];
			p_node_1 = ptcl_num_node[m + 1];
			p_node_2 = ptcl_num_node[m + WidthNum];
			p_node_3 = ptcl_num_node[m + WidthNum + 1];
			p_node_4 = ptcl_num_node[m + WidthNum*DepthNum];
			p_node_5 = ptcl_num_node[m + WidthNum*DepthNum + 1];
			p_node_6 = ptcl_num_node[m + WidthNum*DepthNum + WidthNum];
			p_node_7 = ptcl_num_node[m + WidthNum*DepthNum + WidthNum + 1];

			//対象の要素を構成する全節点に、同じ粒子の情報が記憶されていれば、粒子と要素が重なっている
			if ((p_node_0 == p_node_1) && (p_node_0 == p_node_2) && (p_node_0 == p_node_3) && (p_node_0 == p_node_4) && (p_node_0 == p_node_5) && (p_node_0 == p_node_6) && (p_node_0 == p_node_7)){

				//要素の番号を取得し、粒子情報と粒子の誘電率を記憶する
				ElePerm[m] = spElePerm * ParticleRelElePerm;
				Conductivity[m] = conductivity;
				ptcl_num_ele[m] = ptcl_num_node[m];
				P_C[p].ele.push_back(m);
			}
		}
	}

	/*
	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open("check.fld", ios_base::out | ios_base::trunc);
	fout_charge2.open("check.dat", ios_base::out | ios_base::trunc);

	fout_charge1
	<< "# AVS field file\n" << "ndim = 3\n"
	<< "dim1 = " << WidthNum << '\n'
	<< "dim2 = " << DepthNum << '\n'
	<< "dim3 = " << HeightNum << '\n'
	<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
	<< "field = uniform\n" << '\n'
	<< "variable 1 file=" << "check.dat" << ' ' << "filetype=ascii\n";

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){
				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_charge2 << ptcl_num_ele[m] << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();
	*/

}

void Electric_Calc::Particle_count_charge(){

	cout << "粒子の最外位置に存在する節点を特定します" << endl;

	int n_0, e_1, e_2, e_3, e_4, e_5, e_6;

	//粒子を構成する外周節点にのみ電荷を配置する
	for (int p = 0; p < ptclNum; p++){
		int ele = (int)P_C[p].ele.size();

		for (int d = 0; d < ele; d++){
			int m = P_C[p].ele[d];
			int i = (m % (WidthNum*DepthNum)) % WidthNum;
			int j = (m % (WidthNum*DepthNum)) / WidthNum;
			int k = m / (WidthNum*DepthNum);

			n_0 = ptcl_num_ele[m];
			e_5 = m + WidthNum*DepthNum;
			e_6 = m - WidthNum*DepthNum;

			//注目要素のX・Y方向に隣接する要素情報を取得
			//X方向が0の要素
			if (i == 0){
				//Y方向が0の要素
				if (j == 0){
					e_1 = m + 1;
					e_2 = m + WidthNum - 1;
					e_3 = m + WidthNum;
					e_4 = m + WidthNum*DepthNum - WidthNum;
					//Y方向がMAXの要素
				}
				else if (j == DepthNum - 1){
					e_1 = m + 1;
					e_2 = m + WidthNum - 1;
					e_3 = m - WidthNum*DepthNum + WidthNum;
					e_4 = m - WidthNum;
					//Y方向が端要素以外
				}
				else{
					e_1 = m + 1;
					e_2 = m + WidthNum - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum;
				}
				//X方向がMAXの要素
			}
			else if (i == WidthNum - 1){
				//Y方向が0の要素
				if (j == 0){
					e_1 = m + 1 - WidthNum;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum + WidthNum*DepthNum;
					//Y方向がMAXの要素
				}
				else if (j == DepthNum - 1){
					e_1 = m - WidthNum + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum - WidthNum*DepthNum;
					e_4 = m - WidthNum;
					//Y方向が端要素以外
				}
				else{
					e_1 = m + 1 - WidthNum;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum;
				}
				//X方向が端要素以外
			}
			else{
				//Y方向が0の要素
				if (j == 0){
					e_1 = m + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m + WidthNum*DepthNum - WidthNum;
					//Y方向がMAX
				}
				else if (j == DepthNum - 1){
					e_1 = m + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum - WidthNum*DepthNum;
					e_4 = m - WidthNum;
					//Y方向が端点以外
				}
				else{
					e_1 = m + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum;
				}
			}

			//注目要素と周囲要素に同じ粒子番号が設定されていれば要素は隣接している。よって、最外位置ではない 以下は最外位置を求めている
			if (!(n_0 == ptcl_num_ele[e_1])){
				P_C[n_0].contain.push_back(m + 1);
				P_C[n_0].contain.push_back(m + 1 + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1 + WidthNum);
			}
			if (!(n_0 == ptcl_num_ele[e_2])){
				P_C[n_0].contain.push_back(m);
				P_C[n_0].contain.push_back(m + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + WidthNum);
			}
			if (!(n_0 == ptcl_num_ele[e_3])){
				P_C[n_0].contain.push_back(m + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum + 1);
				P_C[n_0].contain.push_back(m + WidthNum + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum + WidthNum*DepthNum + 1);
			}
			if (!(n_0 == ptcl_num_ele[e_4])){
				P_C[n_0].contain.push_back(m);
				P_C[n_0].contain.push_back(m + 1);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1);
			}
			if (!(n_0 == ptcl_num_ele[e_5])){
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + WidthNum + 1);
			}

			//Z方向の下端面の要素
			if (k == 0){

				P_C[n_0].contain.push_back(m);
				P_C[n_0].contain.push_back(m + 1);
				P_C[n_0].contain.push_back(m + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum + 1);

				//Z方向の端面以外の要素
			}
			else{
				if (!(n_0 == ptcl_num_ele[e_6])){
					P_C[n_0].contain.push_back(m);
					P_C[n_0].contain.push_back(m + 1);
					P_C[n_0].contain.push_back(m + WidthNum);
					P_C[n_0].contain.push_back(m + WidthNum + 1);
				}
			}
		}
	}

	//同じ要素を削除する
	for (int i = 0; i<ptclNum; i++){

		//sortして要素の順番を昇順で並び替える
		sort(P_C[i].contain.begin(), P_C[i].contain.end());
		//uniqueをして連続する同じ要素を連結させ、後ろ残ったゴミをeraseで削除する
		P_C[i].contain.erase(unique(P_C[i].contain.begin(), P_C[i].contain.end()), P_C[i].contain.end());
	}

	/*
	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open("check.fld", ios_base::out | ios_base::trunc);
	fout_charge2.open("check.dat", ios_base::out | ios_base::trunc);

	fout_charge1
	<< "# AVS field file\n" << "ndim = 3\n"
	<< "dim1 = " << WidthNum << '\n'
	<< "dim2 = " << DepthNum << '\n'
	<< "dim3 = " << HeightNum << '\n'
	<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
	<< "field = uniform\n" << '\n'
	<< "variable 1 file=" << "check.dat" << ' ' << "filetype=ascii\n";

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				int temp = 0;

				//sortして要素の順番を昇順で並び替える
				for (int n = 0; n < P_C[0].contain.size(); n++)		if (P_C[0].contain[n] == m)		temp = 10;
				fout_charge2 << temp << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();
	*/

}

void Electric_Calc::Particle_charge_Mapping(){

	//粒子の電荷を配置する
	for (int i = 0; i<ptclNum; i++){

		int Size = int(P_C[i].contain.size());

		for (int j = 0; j<Size; j++){

			int cell_num = P_C[i].contain[j];
			Charge[cell_num] = P_C[i].Ptcl_Charge / Size;
		}
	}
}

void Electric_Calc::Particle_clear(){

	ptcl_num_node.clear();
	ptcl_num_ele.clear();
}


void Electric_Calc::Calc_Potential_Laplace(){

	int localCount = 0;
	double P, P_fcc, P_bcc, P_ccf, P_ccb, P_cfc, P_cbc;
	double e_bbb, e_bbf, e_bfb, e_bff, e_fbb, e_fbf, e_ffb, e_fff;
	int m;

	do{
		#pragma omp barrier
		localCount = 0;

		#pragma omp single
		convergentCount = 0;

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k+=2){
			for (int j = 0; j<DepthNum; j+=2){
				for (int i = 0; i<WidthNum; i+=2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];							
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR 
										* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
											+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
											+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)) 
												/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));
											
						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}					
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp critical
		convergentCount += localCount;

		#pragma omp barrier

#pragma omp single
		cout << convergentCount << "\r";

	} while (convergentCount);
}

void Electric_Calc::Calc_Potential_Poisson(){

	int localCount = 0;
	double Q, P, P_fcc, P_bcc, P_ccf, P_ccb, P_cfc, P_cbc;
	double e_bbb, e_bbf, e_bfb, e_bff, e_fbb, e_fbf, e_ffb, e_fff;
	int m;

	do{
#pragma omp barrier
		localCount = 0;

#pragma omp single
		convergentCount = 0;

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//節点番号を求める
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//節点の電位	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp critical
		convergentCount += localCount;

#pragma omp barrier

	} while (convergentCount);
}

void Electric_Calc::Calc_E(){

	//電界を計算する
	if (Type_boundary_X == 0)		Calc_E_X_n();
	else if (Type_boundary_X == 1)	Calc_E_X_l();

	if (Type_boundary_Y == 0)		Calc_E_Y_n();
	else if (Type_boundary_Y == 1)	Calc_E_Y_l();

	Calc_E_Z_n();
}

void Electric_Calc::Calc_E_X_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (i == 0){
					E_x[m] = 0.0;
				}
				else if (i == WidthNum - 1){
					E_x[m] = 0.0;
				}
				else{
					E_x[m] = -(Potential[m + 1] - Potential[m - 1]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Calc_E_X_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (i == 0){
					E_x[m] = -(Potential[m + 1] - Potential[m + WidthNum - 2]) / (ElmSize * 2.0);
				}
				else if (i == WidthNum - 1){
					E_x[m] = -(Potential[m - WidthNum + 2] - Potential[m - 1]) / (ElmSize * 2.0);
				}
				else{
					E_x[m] = -(Potential[m + 1] - Potential[m - 1]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Calc_E_Y_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (j == 0){
					E_y[m] = 0.0;
				}
				else if (j == DepthNum - 1){
					E_y[m] = 0.0;
				}
				else{
					E_y[m] = -(Potential[m + WidthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Calc_E_Y_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (j == 0){
					E_y[m] = -(Potential[m + WidthNum] - Potential[m + (WidthNum - 2)*DepthNum]) / (ElmSize * 2.0);
				}
				else if (j == DepthNum - 1){
					E_y[m] = -(Potential[m - (WidthNum - 2)*DepthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
				}
				else{
					E_y[m] = -(Potential[m + WidthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Calc_E_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (k == 0){
					E_z[m] = 0.0;
				}
				else if (k == HeightNum - 1){
					E_z[m] = 0.0;
				}
				else{
					E_z[m] = -(Potential[m + (WidthNum*DepthNum)] - Potential[m - (WidthNum*DepthNum)]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Ini_P(){

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				//節点番号を求める
				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				Potential[m] = 0.0;
			}
		}
	}
}

void Electric_Calc::Calc_Ohm(){

	int localCount = 0;
	int n_1, n_2, n_3, n_4, n_5, n_6;
	double e_1, e_2, e_3, e_4, e_5, e_6, e_7, e_8;
	double Jxfwd, Jxbck, Jyfwd, Jybck, Jzfwd, Jzbck;		//各方向の電流密度

	//電荷の計算
	for (int p = 0; p < ptclNum; p++){
		int node = (int)P_C[p].node.size();

		for (int i = 0; i < node; i++){

			int m = P_C[p].node[i];
							
				n_1 = m + 1;
				n_2 = m - 1;
				n_3 = m + WidthNum;
				n_4 = m - WidthNum;
				n_5 = m + WidthNum*DepthNum;
				n_6 = m - WidthNum*DepthNum;

				e_1 = Conductivity[m - WidthNum*DepthNum - WidthNum - 1];
				e_2 = Conductivity[m - WidthNum*DepthNum - WidthNum];
				e_3 = Conductivity[m - WidthNum*DepthNum - 1];
				e_4 = Conductivity[m - WidthNum*DepthNum];
				e_5 = Conductivity[m - WidthNum - 1];
				e_6 = Conductivity[m - WidthNum];
				e_7 = Conductivity[m - 1];
				e_8 = Conductivity[m];

				Jxfwd = -1.0 * (e_2 + e_4 + e_6 + e_8) / 4 * (Potential[n_1] - Potential[m]) / ElmSize;
				Jxbck = -1.0 * (e_1 + e_3 + e_5 + e_7) / 4 * (Potential[m] - Potential[n_2]) / ElmSize;
				Jyfwd = -1.0 * (e_3 + e_4 + e_7 + e_8) / 4 * (Potential[n_3] - Potential[m]) / ElmSize;
				Jybck = -1.0 * (e_1 + e_2 + e_5 + e_6) / 4 * (Potential[m] - Potential[n_4]) / ElmSize;
				Jzfwd = -1.0 * (e_5 + e_6 + e_7 + e_8) / 4 * (Potential[n_5] - Potential[m]) / ElmSize;
				Jzbck = -1.0 * (e_1 + e_2 + e_3 + e_4) / 4 * (Potential[m] - Potential[n_6]) / ElmSize;

				Charge[m] += (-Jxfwd + Jxbck - Jyfwd + Jybck - Jzfwd + Jzbck) * ElmSize * ElmSize * delta_T;
		}			
	}
}

void Electric_Calc::OutPut_Potential(char *filename_1, char *filename_2){

	cout << "ポテンシャルデータを出力します" << endl;

	ofstream fout_potential1;
	ofstream fout_potential2;

	fout_potential1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_potential2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_potential1
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_potential2 << Potential[m] << '\n';
			}
		}
	}

	fout_potential1.close();
	fout_potential2.close();
}

void Electric_Calc::OutPut_E(char *filename_1, char *filename_2){

	cout << "電界データを出力します" << endl;

	ofstream fout_E1;
	ofstream fout_E2;

	fout_E1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_E2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_E1
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_E2 << sqrt(pow(E_x[m], 2) + pow(E_y[m], 2) + pow(E_z[m], 2)) << '\n';
			}
		}
	}

	fout_E1.close();
	fout_E2.close();
}

void Electric_Calc::OutPut_Charge(char *filename_1, char *filename_2){

	cout << "電荷分布データを出力します" << endl;

	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_charge2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_charge1
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_charge2 << Charge[m] * 10e+12 << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();
}

void Electric_Calc::OutPut_Ptcl(char *filename){

	cout << "粒子データを出力します" << endl;
	ofstream fout;
	fout.open(filename, ios_base::out | ios_base::trunc);
	fout << "# Micro AVS Geom:2.00\n";
	fout << "sphere" << endl;
	fout << "Movement" << endl;
	fout << "color" << endl;
	fout << ptclNum << endl;

	for (int i = 0; i < ptclNum; i++){
		fout << P_C[i].Ptcl_Pos_X << " " << P_C[i].Ptcl_Pos_Y << " " << P_C[i].Ptcl_Pos_Z << " " << P_C[i].Ptcl_Radius << " 1 0 0" << endl;
	}

	fout << "disjoint polygon" << endl;
	fout << "floor" << endl;
	fout << "facet" << endl;
	fout << "color" << endl;
	fout << 1 << endl;
	double gray_color = 0.8;

	fout << "4" << endl;
	fout << "0" << " " << "0" << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
	fout << areaWidth << " " << "0" << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
	fout << areaWidth << " " << areaDepth << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
	fout << "0" << " " << areaDepth << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
}

void Electric_Calc::OutPut_Force_1(char *filename){

	ofstream fout;
	fout.open(filename, ios_base::out | ios_base::trunc);

	int P_x_1 = (int)(P_C[0].Ptcl_Pos_X / ElmSize + 0.5);
	int P_y_1 = (int)(P_C[0].Ptcl_Pos_Y / ElmSize + 0.5);
	int P_z_1 = (int)(P_C[0].Ptcl_Pos_Z / ElmSize + 0.5);
	int m_1 = P_z_1 * WidthNum * DepthNum + P_y_1 * WidthNum + P_x_1;

	double E_x_1 = E_x[m_1 + 1] * E_x[m_1 + 1] + E_y[m_1 + 1] * E_y[m_1 + 1] + E_z[m_1 + 1] * E_z[m_1 + 1];
	double E_xx_1 = E_x[m_1 - 1] * E_x[m_1 - 1] + E_y[m_1 - 1] * E_y[m_1 - 1] + E_z[m_1 - 1] * E_z[m_1 - 1];
	double E_y_1 = E_x[m_1 + WidthNum] * E_x[m_1 + WidthNum] + E_y[m_1 + WidthNum] * E_y[m_1 + WidthNum] + E_z[m_1 + WidthNum] * E_z[m_1 + WidthNum];
	double E_yy_1 = E_x[m_1 - WidthNum] * E_x[m_1 - WidthNum] + E_y[m_1 - WidthNum] * E_y[m_1 - WidthNum] + E_z[m_1 - WidthNum] * E_z[m_1 - WidthNum];
	double E_z_1 = E_x[m_1 + WidthNum * DepthNum] * E_x[m_1 + WidthNum * DepthNum] + E_y[m_1 + WidthNum * DepthNum] * E_y[m_1 + WidthNum * DepthNum] + E_z[m_1 + WidthNum * DepthNum] * E_z[m_1 + WidthNum * DepthNum];
	double E_zz_1 = E_x[m_1 - WidthNum * DepthNum] * E_x[m_1 - WidthNum * DepthNum] + E_y[m_1 - WidthNum * DepthNum] * E_y[m_1 - WidthNum * DepthNum] + E_z[m_1 - WidthNum * DepthNum] * E_z[m_1 - WidthNum * DepthNum];

	double temp = 2 * PI * spElePerm * AirRelElePerm * (ParticleRelElePerm - AirRelElePerm) / (ParticleRelElePerm - 2 * AirRelElePerm);

	fout << "粒子1" << endl;
	fout << "時間" << "," << "外力_X" << "," << "外力_Y" << "," << "外力_Z" << endl;
	fout << "0" << "," << P_C[0].Ptcl_Charge * E_x[m_1] << "," << P_C[0].Ptcl_Charge * E_y[m_1] << "," << P_C[0].Ptcl_Charge * E_z[m_1] << endl;

	fout << "粒子1" << endl;
	fout << "時間" << "," << "外力_X" << "," << "外力_Y" << "," << "外力_Z" << endl;
	fout << "0" << "," << temp * (E_x_1 - E_xx_1) / (2.0*ElmSize) * pow(P_C[0].Ptcl_Radius, 3) << "," << temp * (E_y_1 - E_yy_1) / (2.0*ElmSize) * pow(P_C[0].Ptcl_Radius, 3) << "," << temp * (E_z_1 - E_zz_1) / (2.0*ElmSize) * pow(P_C[0].Ptcl_Radius, 3) << endl;

	fout << "以下，ループ" << endl;
	fout << "粒子1" << endl;
	fout << "時間" << "," << "外力_X" << "," << "外力_Y" << "," << "外力_Z" << endl;

	fout.close();
}

void Electric_Calc::OutPut_Force_2(char *filename, int step){

	ofstream fout;
	fout.open(filename, ios_base::app);

	double Force_1_x = 0, Force_1_y = 0, Force_1_z = 0;

	for (int i = 0; i < P_C[0].node.size(); i++){

		int m = P_C[0].node[i];
		Force_1_x += Charge[m] * E_x[m];
		Force_1_y += Charge[m] * E_y[m];
		Force_1_z += Charge[m] * E_z[m];
	}

	fout << step*delta_T << "," << Force_1_x << "," << Force_1_y << "," << Force_1_z << endl;

	fout.close();
}

void Electric_Calc::OutPut_EleFieldBinary(char *filename){

	cout << "挙動計算用の電界データを出力します" << endl;

	ofstream ofs;
	ofs.open(filename, ios_base::out | ios_base::binary | ios_base::trunc);
	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&Potential[m], sizeof(double));
			}
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&E_x[m], sizeof(double));
			}
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&E_y[m], sizeof(double));
			}
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&E_z[m], sizeof(double));
			}
		}
	}
}

void Electric_Calc::GenerateMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat){

	ofstream movie_1;
	ofstream movie_2;
	ofstream movie_3;
	ofstream movie_4;
	ofstream movie_5;
	ofstream movie_6;

	movie_1.open(filename_P_fld, ios_base::out, ios_base::trunc);
	movie_2.open(filename_E_fld, ios_base::out, ios_base::trunc);
	movie_3.open(filename_C_fld, ios_base::out, ios_base::trunc);
	movie_4.open(filename_P_dat, ios_base::out, ios_base::trunc);
	movie_5.open(filename_E_dat, ios_base::out, ios_base::trunc);
	movie_6.open(filename_C_dat, ios_base::out, ios_base::trunc);

	movie_1 << "# AVS field file\n" << "#\n"
		<< "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n"
		<< "veclen = 1\n"
		<< "data = double\n"
		<< "field = uniform\n" << '\n';

	movie_2 << "# AVS field file\n" << "#\n"
		<< "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n"
		<< "veclen = 1\n"
		<< "data = double\n"
		<< "field = uniform\n" << '\n';

	movie_3 << "# AVS field file\n" << "#\n"
		<< "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n"
		<< "veclen = 1\n"
		<< "data = double\n"
		<< "field = uniform\n" << '\n';

	movie_1.close();
	movie_2.close();
	movie_3.close();
	movie_4.close();
	movie_5.close();
	movie_6.close();
}

void Electric_Calc::WriteMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat, int T_step, int Count){

	ofstream movie_1;
	ofstream movie_2;
	ofstream movie_3;
	ofstream movie_4;
	ofstream movie_5;
	ofstream movie_6;

	movie_1.open(filename_P_fld, ios_base::app);
	movie_2.open(filename_E_fld, ios_base::app);
	movie_3.open(filename_C_fld, ios_base::app);
	movie_4.open(filename_P_dat, ios_base::app);
	movie_5.open(filename_E_dat, ios_base::app);
	movie_6.open(filename_C_dat, ios_base::app);

	int STEP;
	STEP = Count*(cellNum_2D + 1);

	int P_y = (int)(P_C[0].Ptcl_Pos_Y / ElmSize + 0.5);

	movie_1 << "time file=" << filename_P_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP << ' ' << "close=1\n"
		<< "variable 1 file=" << filename_P_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP + 1 << ' ' << "close=1\n" << "EOT\n";

	movie_2 << "time file=" << filename_E_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP << ' ' << "close=1\n"
		<< "variable 1 file=" << filename_E_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP + 1 << ' ' << "close=1\n" << "EOT\n";

	movie_3 << "time file=" << filename_C_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP << ' ' << "close=1\n"
		<< "variable 1 file=" << filename_C_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP + 1 << ' ' << "close=1\n" << "EOT\n";

	movie_4 << "time value=" << T_step * delta_T << '\n';
	movie_5 << "time value=" << T_step * delta_T << '\n';
	movie_6 << "time value=" << T_step * delta_T << '\n';

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum * DepthNum + P_y * WidthNum + i;
			movie_4 << Potential[m] << '\n';
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum * DepthNum + P_y * WidthNum + i;
			movie_5 << sqrt(pow(E_x[m], 2) + pow(E_y[m], 2) + pow(E_z[m], 2)) << '\n';

		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum * DepthNum + P_y * WidthNum + i;
			movie_6 << Charge[m] * 10e+12 << '\n';

		}
	}

	movie_1.close();
	movie_2.close();
	movie_3.close();
	movie_4.close();
	movie_5.close();
	movie_6.close();
}

void Electric_Calc::OutPut_Condition(char *filename){

	cout << "計算条件を出力します" << endl;

	ofstream fout_E;
	fout_E.open(filename, ios_base::out | ios_base::trunc);

	fout_E << "mesh_size =  " << ElmSize << " [m] " << endl;
	fout_E << "x_direction =  " << areaWidth << " [m] " << endl;
	fout_E << "y_direction =  " << areaDepth << " [m] " << endl;
	fout_E << "z_direction =  " << areaHeight << " [m] " << endl;
	fout_E << "x_mesh =  " << WidthNum << endl;
	fout_E << "y_mesh =  " << DepthNum << endl;
	fout_E << "z_mesh =  " << HeightNum << endl;
	fout_E << "Volt =  " << Volt << " [V] " << endl;

	fout_E.close();
}