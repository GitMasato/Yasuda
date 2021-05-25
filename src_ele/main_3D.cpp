#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <string>
#include <math.h>
#include <time.h>
#include <omp.h>

using namespace std;
#include "Ele_Calc.h"
#include "Profiler.h"

//main関数部
int main(){

	std::ios::sync_with_stdio(false);

	//入力データの名前		
	char particlefilename[64] = "1ptcl.csv";

	//出力データの名前
	char Ptclfilename[64] = "Ptcl.mgf";							//読み込み用の粒子データ
	char Forcefilename[64] = "Force.csv";						//粒子に加わる外力を書き込むデータ	
	char Boundaryfilename[64] = "Boundary.csv";					//電極上の電位データ	
	char E_Field_1[64] = "1FieldData_3D.bin";					//挙動計算用の電界データ
	char memo[64] = "memo.txt";

	char Chargefilename_1[64] = "Charge.fld";					//電荷データ
	char Chargefilename_2[64] = "Charge.dat";					//電荷データ
	char Laplace_Potentialfilename_1[64] = "Potential.fld";		//ラプラス後の電位データ
	char Laplace_Potentialfilename_2[64] = "Potential.dat";		//ラプラス後の電位データ
	char Laplace_E_filename_1[64] = "ElectricField.fld";		//ラプラス後の電界データ
	char Laplace_E_filename_2[64] = "ElectricField.dat";		//ラプラス後の電界データ

	char P_Moviefilename_1[64] = "P_Movie.fld";					//電位データ	動画
	char P_Moviefilename_2[64] = "P_Movie.dat";					//電位データ	動画		
	char E_Moviefilename_1[64] = "E_Movie.fld";					//電界データ	動画
	char E_Moviefilename_2[64] = "E_Movie.dat";					//電界データ	動画	
	char C_Moviefilename_1[64] = "C_Movie.fld";					//電荷データ	動画
	char C_Moviefilename_2[64] = "C_Movie.dat";					//電荷データ	動画

	//変数の初期化等
	Electric_Calc ele;
	int Num = ele.Get_thread();
	bool OMP_calc = ele.Get_OMP_Calc();
	ele.Setting_Boundary_Condition_V();		//境界条件の設定
	ele.Setting_Boundary_Condition_e();		//境界条件の設定
	ele.Setting_Boundary_Condition_1();		//境界条件の設定

	//計算時間を出力するファイルを作成
	Profiler p;
	ofstream fout;
	fout.open("elapsed_time.txt", ios_base::out, ios_base::trunc);

	if (ele.Get_load_ptcl_Poisson()){

		ele.Ptcl_Data(particlefilename);		//粒子データを読み込み	
		ele.OutPut_Ptcl(Ptclfilename);			//粒子位置を出力	
		ele.Particle_Node_Mapping();			//粒子を節点上に配置	
		ele.Particle_Ele_Mapping();				//粒子を要素上に配置	
		//  ele.Particle_count_charge();		//粒子の最外位置に存在する節点を特定	
		//	ele.Particle_charge_Mapping();		//粒子電荷を最外節点上に配置	
		//	ele.OutPut_Charge(Chargefilename_1,Chargefilename_2);	//粒子電荷データを出力
		ele.Particle_clear();					//不要なメモリを解放

		#pragma omp parallel num_threads(Num) if(OMP_calc)
		{
			ele.Calc_Potential_Poisson();
		}
	}
	else{
		#pragma omp parallel num_threads(Num) if(OMP_calc)
		{
			ele.Calc_Potential_Laplace();
		}
	}

	ele.Calc_E();																		//ラプラス・ポアソン計算後の電界計算	
	ele.OutPut_Potential(Laplace_Potentialfilename_1, Laplace_Potentialfilename_2);		//ラプラス・ポアソンの式を解いた後の電位分布を出力
	ele.OutPut_E(Laplace_E_filename_1, Laplace_E_filename_2);							//ラプラス・ポアソンの式を解いた後の電界分布を出力		
	ele.OutPut_EleFieldBinary(E_Field_1);												//挙動計算用の電界データを出力
	ele.OutPut_Condition(memo);															//計算条件等を出力しておく
	if (ele.Get_load_ptcl_Poisson())	ele.OutPut_Force_1(Forcefilename);				//粒子に加わる外力データの作成　ラプラス・ポアソンのデータを用いて外力を計算

	//時間の出力
	fout << "ラプラスorポアソンの式を解き終わるまで : " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
	cout << "ラプラスorポアソンの式を解き終わるまで : " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;

	if (ele.Get_unsteady_calc()){

		//電位分布を初期化
		ele.Ini_P();

		//動画ファイルの生成
		cout << "GenerateMovieFile" << endl;
		ele.GenerateMovieFile(P_Moviefilename_1, P_Moviefilename_2, E_Moviefilename_1, E_Moviefilename_2, C_Moviefilename_1, C_Moviefilename_2);

		//動画の書き込み					
		ele.WriteMovieFile(P_Moviefilename_1, P_Moviefilename_2, E_Moviefilename_1, E_Moviefilename_2, C_Moviefilename_1, C_Moviefilename_2, 0, 0);

		int Counter_1 = 1, Counter_2 = 0;
		int WriteAVS = ele.Get_AVS();

		//ループ開始
		cout << "Calculating" << endl;
		for (int step = 0; step < ele.Return_step(); step++){

			ele.Setting_Boundary_Condition_2(step, Boundaryfilename);		//境界条件の修正

			if (ele.Get_load_ptcl_Poisson()){

				#pragma omp parallel num_threads(Num) if(OMP_calc)
				{
					ele.Calc_Potential_Poisson();
				}
				ele.Calc_E();								//電界計算			
				ele.OutPut_Force_2(Forcefilename, step);	//外力の計算			
				ele.Calc_Ohm();								//電荷保存則の計算
			}
			else{
				#pragma omp parallel num_threads(Num) if(OMP_calc)
				{
					ele.Calc_Potential_Laplace();
				}
				ele.Calc_E();								//電界計算
			}

			//動画の書き込み　時間の出力
			if (WriteAVS == Counter_2){

				ele.WriteMovieFile(P_Moviefilename_1, P_Moviefilename_2, E_Moviefilename_1, E_Moviefilename_2, C_Moviefilename_1, C_Moviefilename_2, step, Counter_1);
				fout << " step = " << step << " / " << ele.Return_step() << " erapedtime = : " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
				cout << " step = " << step << " / " << ele.Return_step() << " erapedtime = : " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
				Counter_1++;
				Counter_2 = 0;
			}
			Counter_2++;			
		}
	}	
}