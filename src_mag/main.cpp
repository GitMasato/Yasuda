#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <string>
#include <cmath>
#include <omp.h>

using namespace std;
#include "Mag_Calc.h"
#include "Profiler.h"

//main関数部
int main(){

	std::ios::sync_with_stdio(false);
	Profiler p;

	//exe等の実行ファイルと同じフォルダにOUTPUTフォルダを作成しておく
	//出力ファイルの名前
	char B_Field_0[64] = "1Mag_FieldData_3D_0.bin";				//ポアソン後の粒子挙動計算用データ
	char B_Field_1[64] = "1Mag_FieldData_3D_1.bin";				//ポアソン後の粒子挙動計算用データ
	char B_Field_2[64] = "1Mag_FieldData_3D_2.bin";				//ポアソン後の粒子挙動計算用データ
	char B_Field_3[64] = "1Mag_FieldData_3D_3.bin";				//ポアソン後の粒子挙動計算用データ

	char Cartesian_0_a[64] = "Cartesian_0.fld";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ
	char Cartesian_0_b[64] = "Cartesian_0.dat";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ
	char Cartesian_1_a[64] = "Cartesian_1.fld";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ
	char Cartesian_1_b[64] = "Cartesian_1.dat";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ
	char Cartesian_2_a[64] = "Cartesian_2.fld";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ
	char Cartesian_2_b[64] = "Cartesian_2.dat";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ
	char Cartesian_3_a[64] = "Cartesian_3.fld";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ
	char Cartesian_3_b[64] = "Cartesian_3.dat";						//ポアソン後の粒子挙動計算用データをチェックするためのデータ

	char Poisson_Potentialfilename_1[64] = "P_Poisson.fld";		//ポアソン後のベクトルポテンシャルデータ
	char Poisson_Potentialfilename_2[64] = "P_Poisson.dat";		//ポアソン後のベクトルポテンシャルデータ
	char Poisson_B_filename_1[64] = "B_Poisson_1.fld";			//ポアソン後の磁束密度データ
	char Poisson_B_filename_2[64] = "B_Poisson_1.dat";			//ポアソン後の磁束密度データ
	char Poisson_B_filename_3[64] = "B_Poisson_2.fld";			//ポアソン後の磁束密度データ
	char Poisson_B_filename_4[64] = "B_Poisson_2.dat";			//ポアソン後の磁束密度データ
	char Poisson_Br[64] = "Br.csv";                             //コイルの半分の位置の磁束密度データ
	char Poisson_Bz[64] = "Bz.csv";								//中心線上の磁束密度データ
	
	char Coil_filename_1[64] = "Coil_Position.fld";				//コイル位置データ
	char Coil_filename_2[64] = "Coil_Position.dat";				//コイルデータ
	
	char Peidong_Bz_0[64] = "Bz_dis_0.csv";						//一時的なデータ
	char Peidong_Bz_1[64] = "Bz_dis_1.csv";						//一時的なデータ
	char Peidong_Bz_2[64] = "Bz_dis_2.csv";						//一時的なデータ
	char Peidong_Bz_3[64] = "Bz_dis_3.csv";						//一時的なデータ
	char memo[64] = "memo.txt";
	
	
	
	//変数の初期化等
	Magnetic_Calc mag;
	int Num = mag.Get_thread();
	bool OMP_calc = mag.Get_OMP_Calc();

	//コイル1巻の位置など計算
	mag.Calc_CoilPosition();

	//境界条件の設定
	#pragma omp parallel num_threads(Num) if(OMP_calc)
	{
		mag.Setting_Boundary_Condition();
	}

	//コイル位置を出力
	mag.OutPut_Coil(Coil_filename_1, Coil_filename_2);	
	
	//ポアソンの式
	#pragma omp parallel num_threads(Num) if(OMP_calc)
	{
		mag.Calc_Potential_Poisson();
	}

	//ポアソンの式を解いた後のベクトルポテンシャル分布を出力
	mag.OutPut_Potential(Poisson_Potentialfilename_1, Poisson_Potentialfilename_2);
	
	//磁束密度計算
	mag.Calc_B();				

	//ポアソンの式を解いた後の磁束密度分布を出力
	mag.OutPut_B_1(Poisson_B_filename_1, Poisson_B_filename_2);
	mag.OutPut_B_2(Poisson_B_filename_3, Poisson_B_filename_4);
	mag.OutPut_B(Poisson_Br,Poisson_Bz);

	//磁束密度データを円筒座標系からデカルト座標へ変換に変換
	mag.Cal_B_Cylindrical();
	mag.Cal_Temp(Peidong_Bz_0, Peidong_Bz_1, Peidong_Bz_2, Peidong_Bz_3);

	//デカルト座標系の磁束密度が正しく出力されているか確認
	mag.OutPut_B_Cylindrical(Cartesian_0_a, Cartesian_0_b, 
		Cartesian_1_a, Cartesian_1_b, 
		Cartesian_2_a, Cartesian_2_b, 
		Cartesian_3_a, Cartesian_3_b);

	//挙動計算用の電界データを出力
	mag.OutPut_MagFieldBinary(B_Field_0, B_Field_1, B_Field_2, B_Field_3);

	//計算条件等を出力
	mag.OutPut_Condition(memo);
	
	//計算時間を出力する
	ofstream fout;
	fout.open("2_elapsed_time.txt", ios_base::out | ios_base::trunc);
	fout << "erapsed time = " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
	fout.close();
	
}