#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <string>
#include <cmath>
#include <ctime>
#include <omp.h>

using namespace std;
#include "Ele_Calc.h"
#include "Profiler.h"

//main関数部
int main(){
	
	std::ios::sync_with_stdio(false);
	
	//exe等の実行ファイルと同じフォルダにOUTPUTフォルダを作成しておく
	//出力ファイルの名前
	char E_Field_1[64] = "1FieldData_2D.bin";						//ラプラス後のデータ
	char End_Laplace_Potentialfilename_1[64] = "P_Laplace.fld";		//ラプラス後の電位データ
	char End_Laplace_Potentialfilename_2[64] = "P_Laplace.dat";		//ラプラス後の電位データ
	char End_Laplace_E_filename_1[64] = "E_Laplace.fld";			//ポアソン後の電位データ
	char End_Laplace_E_filename_2[64] = "E_Laplace.dat";			//ポアソン後の電位データ
	char memo[64] = "memo.txt";
		
	//変数の初期化等
	Electric_Calc ele;
	Profiler p;
	int Num = ele.Get_thread();
	bool OMP_calc = ele.Get_OMP_Calc();

	//境界条件の設定
	ele.Setting_Boundary_Condition_V();
	ele.Setting_Boundary_Condition_e();
	ele.Setting_Boundary_Condition();

	//ラプラスの式
	#pragma omp parallel num_threads(Num) if(OMP_calc)
	{
		ele.Calc_Potential_Laplace();
	}

	//ラプラスの式を解いた後の電位分布を出力
	ele.OutPut_Potential(End_Laplace_Potentialfilename_1, End_Laplace_Potentialfilename_2);

	//ラプラス計算後の電界計算
	ele.Calc_E();

	//ラプラスの式を解いた後の電界分布を出力
	ele.OutPut_E(End_Laplace_E_filename_1, End_Laplace_E_filename_2);
	
	//挙動計算用の電界データを出力
	ele.OutPut_EleFieldBinary(E_Field_1);

	//計算条件等を出力しておく
	ele.OutPut_Condition(memo);

	//計算時間を出力する
	ofstream fout;
	fout.open("2_elapsed_time.txt", ios_base::out | ios_base::trunc);
	fout << "erapsed time = " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
	fout.close();

}