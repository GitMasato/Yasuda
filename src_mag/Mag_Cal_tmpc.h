#ifndef MAG_CALC_H		//インクルードガード
#define MAG_CALC_H

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

struct Electrode_info{
	int pos_R;					//位置節点 
	int pos_Z;					//位置節点 
};

struct OutPut_Coil_info{

	int Coil_Pos_X;				// コイル位置 X 方向
	int Coil_Pos_Y;				// コイル位置 Y 方向
	int Coil_Pos_Z;				// コイル位置 Z 方向
	double coil_pos_x;			// コイル位置　X 方向
	double coil_pos_y;			// コイル位置　Y 方向
	double coil_pos_z;			// コイル位置　Z 方向
	double incline_XZ;			// コイルの傾き角度
	double incline_XY;			// コイルの傾き角度
	double SIN_XY;				// コイルの傾き角度
	double COS_XY;				// コイルの傾き角度
	double SIN_XZ;				// コイルの傾き角度
	double COS_XZ;				// コイルの傾き角度

};

class Magnetic_Calc{

private:

	//基本パラメータ
	const double PI;			//円周率
	int wire_turn;				//コイル巻き数
	int turn_number_height;		//コイル長さ方向巻き数
	int nodeNum;				//節点の総数	
	int WidthNum;				// r 方向節点数
	int HeightNum;				// z 方向節点数
	int convergentCount;		//収束していない節点数
	int thread;					//並列化のスレッド数	
	bool OMP_CALC;				//並列化するか
	double ElmSize;				//要素長さ
	double areaWidth;			//計算領域R方向長さ
	double areaHeight;			//計算領域Z方向長さ
	double permeability;		//真空の透磁率
	double wire_diameter;		//導線直径     m
	double wire_totallength;	//導線長さ     m
	double coil_length;			//コイル長さ 
	double coil_diameter;		//コイルを巻く円筒の直径　 m
	double endpoint;			//計算領域下端からコイル先端までの距離　 m
	double distance_coils;		//コイル間の距離
	double I;					//印加電流
	double W;					//消費電力    W
	double p;					//電気抵抗率　Ωm
	double R;					//コイル全体の抵抗 Ω
	double TFC;					//収束判定
	double SOR;					//加速係数		

	//配列データ
	vector<double> Potential;		//節点のベクトルポテンシャル
	vector<double> CurrentDensity;	//節点の電流密度
	vector<double> B_r;				//節点の磁束密度
	vector<double> B_z;				//節点の磁束密度
	vector<double> coefficient_1;	//ベクトルポテンシャル計算用の係数
	vector<double> coefficient_2;	//ベクトルポテンシャル計算用の係数
	vector<double> coefficient_3;	//ベクトルポテンシャル計算用の係数
	vector<double> coefficient_4;	//ベクトルポテンシャル計算用の係数
	bool * CalcFlag;				//節点をポテンシャル計算するかのフラグ
	vector<Electrode_info> electrode_1;
	vector<Electrode_info> electrode_2;
	vector<OutPut_Coil_info> output_coil;

	//粒子挙動計算用の出力パラメータ
	int nodeNum_output;			//節点の総数	
	int WidthNum_output;		// X 方向節点数
	int DepthNum_output;		// Y 方向節点数
	int HeightNum_output;		// Z 方向節点数
	double ElmSize_output;		//要素長さ
	double areaWidth_output;	//計算領域 X 方向長さ　　円筒座標系の　コイル下端部の中心位置から
	double areaDepth_output;	//計算領域 Y 方向長さ　
	double areaHeight_output;	//計算領域 Z 方向長さ	
	vector<double> Output_B_x_0;	//節点の磁束密度　円筒座標
	vector<double> Output_B_y_0;	//節点の磁束密度　円筒座標
	vector<double> Output_B_z_0;	//節点の磁束密度　円筒座標
	vector<double> Output_B_x_1;	//節点の磁束密度　円筒座標
	vector<double> Output_B_y_1;	//節点の磁束密度　円筒座標
	vector<double> Output_B_z_1;	//節点の磁束密度　円筒座標
	vector<double> Output_B_x_2;	//節点の磁束密度　円筒座標
	vector<double> Output_B_y_2;	//節点の磁束密度　円筒座標
	vector<double> Output_B_z_2;	//節点の磁束密度　円筒座標
	vector<double> Output_B_x_3;	//節点の磁束密度　円筒座標
	vector<double> Output_B_y_3;	//節点の磁束密度　円筒座標
	vector<double> Output_B_z_3;	//節点の磁束密度　円筒座標
	

public:

	Magnetic_Calc();							//コンストラクタ
	~Magnetic_Calc();							//デストラクタ
		
	int Get_thread(){ return thread; }			//スレッド数を返す
	bool Get_OMP_Calc(){ return OMP_CALC; }		//並列化するか
	void Calc_CoilPosition();					//コイルの位置特定
	void Setting_Boundary_Condition();			//境界条件の設定
	void Calc_Potential_Poisson();				//ポアソンを解く　ベクトルポテンシャルを求める
	void Calc_B();								//磁束密度計算
	void Cal_B_Cylindrical();				    //磁束密度データを円筒座標系からデカルト座標へ変換に変換
	void Cal_Temp(char *filename_0, char *filename_1, char *filename_2, char *filename_3);

	void OutPut_Coil(char *filename_1, char *filename_2);			//コイル位置を出力する
	void OutPut_Potential(char *filename_1, char *filename_2);		//ベクトルポテンシャルデータを出力する
	void OutPut_B_1(char *filename_1, char *filename_2);			//磁界データを出力する
	void OutPut_B_2(char *filename_1, char *filename_2);			//磁界データを出力する
	void OutPut_B(char *filename_1, char *filename_2);				//磁界データを出力する

	void OutPut_B_Cylindrical(char *filename_0_a, char *filename_0_b,
		char *filename_1_a, char *filename_1_b, 
		char *filename_2_a, char *filename_2_b,
		char *filename_3_a, char *filename_3_b);	//デカルト座標系の磁束密度が正しく出力されているか確認	
	
	
	void OutPut_MagFieldBinary(char *filename_0, char *filename_1, char *filename_2, char *filename_3);	//挙動計算用の磁界データを出力する							
	void OutPut_Condition(char *filename);							//計算条件を出力する

};
#endif