#ifndef ELE_CALC_H		//インクルードガード
#define ELE_CALC_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <cstdlib>

using namespace std;

struct Ptcl_info{

	vector<int> contain;
	vector<int> node;						
	vector<int> ele;		
	double Ptcl_Pos_X;		//粒子の位置
	double Ptcl_Pos_Y;		//粒子の位置
	double Ptcl_Pos_Z;		//粒子の位置
	double Ptcl_Radius;		//粒子の半径
	double Ptcl_Charge;		//粒子の帯電量
	double sp_G;			//比重	
};

struct Cylinder_info{

	int pos_X;						//位置節点 
	int pos_Y;						//位置節点 
	int pos_Z;						//位置節点 
	int radius;
	int length;
	int diele_radius;
	double Pos_X;					//一端点の座標
	double Pos_Y;					//一端点の座標
	double Pos_Z;					//一端点の座標
	double Radius;					//半径
	double Length;					//長さ
	double Dielectric_Radius;		//誘電体
	vector<int> cyl_mesh;	    	//位置節点 
};

struct Hollow_Cylinder_info{

	int pos_X;						//位置節点 
	int pos_Y;						//位置節点 
	int pos_Z;						//位置節点 
	int inside_radius;
	int outside_radius;
	int length;
	double Pos_X;					//端点の座標
	double Pos_Y;					//端点の座標
	double Pos_Z;					//端点の座標
	double Inside_Radius;			//半径
	double Outside_Radius;			//半径
	double Length;					//長さ
	vector<int> hollow_cyl_mesh;	//位置節点 
};

struct Wall_info{

	int pos_X[4];					//位置節点 
	int pos_Y[4];					//位置節点 
	int pos_Z[4];					//位置節点 
	int length;						//長さ
	double Pos_X[4];				//端点の座標
	double Pos_Y[4];				//端点の座標
	double Pos_Z[4];				//端点の座標
	double Length;					//長さ
	vector<int> wall_mesh;			//位置節点 
};

class Electric_Calc{

private:
		
	//基本パラメータ
	int Type_boundary_X;						//境界条件のパターン
	int Type_boundary_Y;						//境界条件のパターン
	int cellNum_2D;								//節点の総数
	int cellNum_3D;								//節点の総数	
	int WidthNum;								//x方向節点数
	int DepthNum;								//y方向節点数
	int HeightNum;								//z方向節点数	
	int convergentCount;						//収束していない節点数
	int thread;									//並列化のスレッド数
	bool OMP_CALC;								//並列化するか	
	bool unsteady_calc;							//非定常の計算か
	bool load_ptcl_Poisson;						//粒子データを読み込むか　ポアソンか
	const double PI;							//円周率

	double ElmSize;								//要素長さ
	double areaWidth;							//計算領域X方向長さ
	double areaDepth;							//計算領域Y方向長さ
	double areaHeight;							//計算領域Z方向長さ
	double spElePerm;							//真空の誘電率
	double AirRelElePerm;						//空気比誘電率
	double Volt;								//印加電圧	
	double TFC;									//収束判定
	double SOR;									//加速係数

	vector<int> bcc;					        //計算する節点	(電位特定用)
	vector<int> fcc;						    //計算する節点
	vector<int> cbc;					        //計算する節点
	vector<int> cfc;						    //計算する節点
	vector<int> ccb;							//計算する節点
	vector<int> ccf;							//計算する節点
	vector<int> bbb;					        //計算する節点	(誘電率特定用)
	vector<int> fbb;						    //計算する節点
	vector<int> bfb;					        //計算する節点
	vector<int> ffb;						    //計算する節点
	vector<int> bbf;					        //計算する節点	
	vector<int> fbf;						    //計算する節点
	vector<int> bff;					        //計算する節点
	vector<int> fff;						    //計算する節点
	vector<double> Potential;					//節点の電位
	vector<double> E_x;							//節点の電界
	vector<double> E_y;							//節点の電界
	vector<double> E_z;							//節点の電界
	vector<double> ElePerm;						//要素の誘電率
	bool * CalcFlag;							//節点をポテンシャル計算するかのフラグ

	//非定常電界計算を行う場合のみ
	int total_step;								//総ステップ
	int phase;									//何相か
	int WriteAVS;								//AVS動画の分割数
	int phase_No;								//プラスの電圧				
	bool phase_change_flag;						//電界を切り替えるか
	double span;								//一時変数
	double delta_T;								//離散化時間
	double total_T;								//計算時間
	double frequency;							//周波数
	double rise_time;							//立ち上がり時間
		
	//粒子データ	ポアソンの式を解く場合のみ	
	int ptclNum;								//粒子数
	double conductivity;						//粒子の導電率
	double ParticleRelElePerm;					//粒子比誘電率	
	vector<Ptcl_info> P_C;						//Ptcl_infoにアクセスするためのインスタンス	
	vector<int> ptcl_num_node;					//節点に配置される粒子の番号			
	vector<int> ptcl_num_ele;					//要素に配置される粒子の番号	
	vector<double> Conductivity;				//要素の導電率	
	vector<double> Charge;						//節点の電荷
	
	//円筒を入れる場合
	int cylinder_Xdir_Num;						
	int cylinder_Ydir_Num;						
	vector<Cylinder_info> cylinder_Xdir;
	vector<Cylinder_info> cylinder_Ydir;
	double CylinderElePerm;

	//中空円筒を入れる場合
	int hollow_cylNum;
	vector<Hollow_Cylinder_info> hollow_cylinder;
	double hollow_cylinderElePerm;

	//直方体を入れる場合
	int wallNum;
	vector<Wall_info> wall;
	double wallElePerm;
		
public:

	Electric_Calc();												//コンストラクタ
	~Electric_Calc();												//デストラクタ

	int Get_AVS(){ return WriteAVS; }								//スレッド数を返す
	int Get_thread(){ return thread; }								//スレッド数を返す
	bool Get_OMP_Calc(){ return OMP_CALC; }							//並列化するか
	bool Get_unsteady_calc(){ return unsteady_calc; }				//非定常の計算か
	bool Get_load_ptcl_Poisson(){ return load_ptcl_Poisson; }		//粒子を考慮するか・ポアソンの式か
	int Return_step(){ return total_step; }							//計算ステップ数を返す関数

	void Setting_Boundary_Condition_V();							//境界条件の設定
	void Setting_Boundary_Condition_V_X_n();						//境界条件の設定
	void Setting_Boundary_Condition_V_X_l();						//境界条件の設定
	void Setting_Boundary_Condition_V_Y_n();						//境界条件の設定
	void Setting_Boundary_Condition_V_Y_l();						//境界条件の設定
	void Setting_Boundary_Condition_V_Z_n();						//境界条件の設定

	void Setting_Boundary_Condition_e();							//境界条件の設定
	void Setting_Boundary_Condition_e_X_n_Y_n_Z_n();				//境界条件の設定

	void Setting_Boundary_Condition_e_X_n_corner();					//境界条件の設定
	void Setting_Boundary_Condition_e_Y_n_corner();					//境界条件の設定
	void Setting_Boundary_Condition_e_Z_n_corner();					//境界条件の設定
	void Setting_Boundary_Condition_e_X_n_Y_n_corner();				//境界条件の設定
	void Setting_Boundary_Condition_e_X_n_Z_n_corner();				//境界条件の設定
	void Setting_Boundary_Condition_e_Y_n_Z_n_corner();				//境界条件の設定
	void Setting_Boundary_Condition_e_X_n_Y_n_Z_n_corner();			//境界条件の設定
	
	void Setting_Boundary_Condition_1();							//境界条件の設定
	void Setting_Boundary_Condition_2(int step, char *filename);	//境界条件の設定
	void Ptcl_Data(char *filename);									//粒子データの読み込み		
	void Particle_Node_Mapping();									//粒子を節点上に設置
	void Particle_Ele_Mapping();									//粒子を要素上に設置
	void Particle_count_charge();									//粒子ごとに、電荷を配置する節点数を数える
	void Particle_charge_Mapping();									//電荷を配置する
	void Particle_clear();											//不要なメモリを解放する
	void Calc_Potential_Laplace();									//ラプラスを解く		
	void Calc_Potential_Poisson();									//ポアソンを解く
	void Calc_E();													//電界計算
	void Calc_E_X_n();												//電界計算
	void Calc_E_X_l();												//電界計算
	void Calc_E_Y_n();												//電界計算
	void Calc_E_Y_l();												//電界計算
	void Calc_E_Z_n();												//電界計算

	void Calc_Ohm();												//オームの式計算		
	void Ini_P();													//電位分布を初期化
	void OutPut_Potential(char *filename_1, char *filename_2);		//電位データを出力する
	void OutPut_E(char *filename_1, char *filename_2);				//電界データを出力する
	void OutPut_Charge(char *filename_1, char *filename_2);			//電荷データを出力する
	void OutPut_Ptcl(char *filename);								//粒子位置を出力する
	void OutPut_Force_1(char *filename);							//粒子に加わる外力を出力する
	void OutPut_Force_2(char *filename, int step);					//粒子に加わる外力を出力する
	void OutPut_EleFieldBinary(char *filename);						//挙動計算用の電界データを出力する
	void GenerateMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat);	//動画データを作成		
	void WriteMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat, int T_step, int Count);	//動画データを書き込み
	void OutPut_Condition(char *filename);

	int int_string(const string& str){			//string型をint型へ変換する関数

		int t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	double double_string(const string& str){			//string型をdouble型へ変換する関数

		double t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	list<string> split(string str, string delim){

		list<string> result;
		int cutAt;

		while ((cutAt = int(str.find_first_of(delim))) != str.npos){

			if (cutAt > 0){

				result.push_back(str.substr(0, cutAt));
			}

			str = str.substr(cutAt + 1);

		}

		if (str.length() > 0){

			result.push_back(str);
		}

		return result;
	}
};
#endif