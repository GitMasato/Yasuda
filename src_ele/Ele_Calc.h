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

struct Cylinder_info{

	int pos_X;						//位置節点 
	int pos_Y;						//位置節点 
	int pos_Z;						//位置節点 
	int radius;
	int diele_radius;
	double Pos_X;					//一端点の座標
	double Pos_Y;					//一端点の座標
	double Pos_Z;					//一端点の座標
	double Radius;					//半径
	double Dielectric_Radius;		//誘電体
};

struct Wall_info{

	int pos_X[4];					//位置節点 
	int pos_Y[4];					//位置節点 
	int pos_Z[4];					//位置節点 
	double Pos_X[4];				//端点の座標
	double Pos_Y[4];				//端点の座標
	double Pos_Z[4];				//端点の座標
};

class Electric_Calc{

private:

	//基本パラメータ
	const double PI;							//円周率
	int Type_boundary_X;						//境界条件のパターン
	int cellNum;								//節点の総数	
	int WidthNum;								//x方向節点数
	int HeightNum;								//z方向節点数
	int convergentCount;						//収束していない節点数
	int thread;									//並列化のスレッド数	
	bool OMP_CALC;								//並列化するか
	double ElmSize;								//要素長さ
	double areaWidth;							//計算領域X方向長さ
	double areaHeight;							//計算領域Z方向長さ
	double spElePerm;							//真空の誘電率．
	double AirRelElePerm;						//空気比誘電率
	double Volt;								//印加電圧
	double TFC;									//収束判定
	double SOR;									//加速係数

	vector<int> bcc;					        //計算する節点の左節点(電位特定用)
	vector<int> fcc;						    //計算する節点の右節点
	vector<int> ccb;							//計算する節点の上節点
	vector<int> ccf;							//計算する節点の下節点
	vector<int> bcb;							//計算する節点の左下節点(誘電率特定用)
	vector<int> bcf;							//計算する節点の右下節点
	vector<int> fcb;							//計算する節点の左上節点
	vector<int> fcf;							//計算する節点の右上節点

	vector<double> ElePerm;						//要素の誘電率
	vector<double> Potential;					//節点の電位
	vector<double> E_x;							//節点の電界
	vector<double> E_z;							//節点の電界
	bool * CalcFlag;							//節点をポテンシャル計算するかのフラグ

	//円を入れる場合
	int cylinder_Num;
	vector<Cylinder_info> cylinder;
	double CylinderElePerm;
	
	//四角を入れる場合
	int wallNum_electrode;
	vector<Wall_info> wall_electrode;
	int wallNum_dielectric;
	vector<Wall_info> wall_dielectric;
	double wallElePerm;

public:

	Electric_Calc();												//コンストラクタ
	~Electric_Calc();												//デストラクタ

	int Get_thread(){ return thread; }								//スレッド数を返す
	bool Get_OMP_Calc(){ return OMP_CALC; }							//並列化するか
	void Setting_Boundary_Condition_V();							//境界条件の設定
	void Setting_Boundary_Condition_V_X_n();						//境界条件の設定
	void Setting_Boundary_Condition_V_X_r();						//境界条件の設定
	void Setting_Boundary_Condition_V_Z_n();						//境界条件の設定

	void Setting_Boundary_Condition_e();							//境界条件の設定
	void Setting_Boundary_Condition_e_X_n_Z_n();					//境界条件の設定
	void Setting_Boundary_Condition_e_X_r_Z_n();					//境界条件の設定

	void Setting_Boundary_Condition();								//境界条件の設定
	void Calc_Potential_Laplace();									//ラプラスを解く
	void Calc_E();													//電界計算
	void Calc_E_X_n();												//電界計算
	void Calc_E_X_r();												//電界計算
	void Calc_E_Z_n();												//電界計算

	void OutPut_Potential(char *filename_1, char *filename_2);		//電位データを出力する
	void OutPut_E(char *filename_1, char *filename_2);				//電界データを出力する
	void OutPut_EleFieldBinary(char *filename);						//挙動計算用の電界データを出力する							
	void OutPut_Condition(char *filename);
	
};
#endif