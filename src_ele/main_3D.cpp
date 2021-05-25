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

//main�֐���
int main(){

	std::ios::sync_with_stdio(false);

	//���̓f�[�^�̖��O		
	char particlefilename[64] = "1ptcl.csv";

	//�o�̓f�[�^�̖��O
	char Ptclfilename[64] = "Ptcl.mgf";							//�ǂݍ��ݗp�̗��q�f�[�^
	char Forcefilename[64] = "Force.csv";						//���q�ɉ����O�͂��������ރf�[�^	
	char Boundaryfilename[64] = "Boundary.csv";					//�d�ɏ�̓d�ʃf�[�^	
	char E_Field_1[64] = "1FieldData_3D.bin";					//�����v�Z�p�̓d�E�f�[�^
	char memo[64] = "memo.txt";

	char Chargefilename_1[64] = "Charge.fld";					//�d�׃f�[�^
	char Chargefilename_2[64] = "Charge.dat";					//�d�׃f�[�^
	char Laplace_Potentialfilename_1[64] = "Potential.fld";		//���v���X��̓d�ʃf�[�^
	char Laplace_Potentialfilename_2[64] = "Potential.dat";		//���v���X��̓d�ʃf�[�^
	char Laplace_E_filename_1[64] = "ElectricField.fld";		//���v���X��̓d�E�f�[�^
	char Laplace_E_filename_2[64] = "ElectricField.dat";		//���v���X��̓d�E�f�[�^

	char P_Moviefilename_1[64] = "P_Movie.fld";					//�d�ʃf�[�^	����
	char P_Moviefilename_2[64] = "P_Movie.dat";					//�d�ʃf�[�^	����		
	char E_Moviefilename_1[64] = "E_Movie.fld";					//�d�E�f�[�^	����
	char E_Moviefilename_2[64] = "E_Movie.dat";					//�d�E�f�[�^	����	
	char C_Moviefilename_1[64] = "C_Movie.fld";					//�d�׃f�[�^	����
	char C_Moviefilename_2[64] = "C_Movie.dat";					//�d�׃f�[�^	����

	//�ϐ��̏�������
	Electric_Calc ele;
	int Num = ele.Get_thread();
	bool OMP_calc = ele.Get_OMP_Calc();
	ele.Setting_Boundary_Condition_V();		//���E�����̐ݒ�
	ele.Setting_Boundary_Condition_e();		//���E�����̐ݒ�
	ele.Setting_Boundary_Condition_1();		//���E�����̐ݒ�

	//�v�Z���Ԃ��o�͂���t�@�C�����쐬
	Profiler p;
	ofstream fout;
	fout.open("elapsed_time.txt", ios_base::out, ios_base::trunc);

	if (ele.Get_load_ptcl_Poisson()){

		ele.Ptcl_Data(particlefilename);		//���q�f�[�^��ǂݍ���	
		ele.OutPut_Ptcl(Ptclfilename);			//���q�ʒu���o��	
		ele.Particle_Node_Mapping();			//���q��ߓ_��ɔz�u	
		ele.Particle_Ele_Mapping();				//���q��v�f��ɔz�u	
		//  ele.Particle_count_charge();		//���q�̍ŊO�ʒu�ɑ��݂���ߓ_�����	
		//	ele.Particle_charge_Mapping();		//���q�d�ׂ��ŊO�ߓ_��ɔz�u	
		//	ele.OutPut_Charge(Chargefilename_1,Chargefilename_2);	//���q�d�׃f�[�^���o��
		ele.Particle_clear();					//�s�v�ȃ����������

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

	ele.Calc_E();																		//���v���X�E�|�A�\���v�Z��̓d�E�v�Z	
	ele.OutPut_Potential(Laplace_Potentialfilename_1, Laplace_Potentialfilename_2);		//���v���X�E�|�A�\���̎�����������̓d�ʕ��z���o��
	ele.OutPut_E(Laplace_E_filename_1, Laplace_E_filename_2);							//���v���X�E�|�A�\���̎�����������̓d�E���z���o��		
	ele.OutPut_EleFieldBinary(E_Field_1);												//�����v�Z�p�̓d�E�f�[�^���o��
	ele.OutPut_Condition(memo);															//�v�Z���������o�͂��Ă���
	if (ele.Get_load_ptcl_Poisson())	ele.OutPut_Force_1(Forcefilename);				//���q�ɉ����O�̓f�[�^�̍쐬�@���v���X�E�|�A�\���̃f�[�^��p���ĊO�͂��v�Z

	//���Ԃ̏o��
	fout << "���v���Xor�|�A�\���̎��������I���܂� : " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
	cout << "���v���Xor�|�A�\���̎��������I���܂� : " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;

	if (ele.Get_unsteady_calc()){

		//�d�ʕ��z��������
		ele.Ini_P();

		//����t�@�C���̐���
		cout << "GenerateMovieFile" << endl;
		ele.GenerateMovieFile(P_Moviefilename_1, P_Moviefilename_2, E_Moviefilename_1, E_Moviefilename_2, C_Moviefilename_1, C_Moviefilename_2);

		//����̏�������					
		ele.WriteMovieFile(P_Moviefilename_1, P_Moviefilename_2, E_Moviefilename_1, E_Moviefilename_2, C_Moviefilename_1, C_Moviefilename_2, 0, 0);

		int Counter_1 = 1, Counter_2 = 0;
		int WriteAVS = ele.Get_AVS();

		//���[�v�J�n
		cout << "Calculating" << endl;
		for (int step = 0; step < ele.Return_step(); step++){

			ele.Setting_Boundary_Condition_2(step, Boundaryfilename);		//���E�����̏C��

			if (ele.Get_load_ptcl_Poisson()){

				#pragma omp parallel num_threads(Num) if(OMP_calc)
				{
					ele.Calc_Potential_Poisson();
				}
				ele.Calc_E();								//�d�E�v�Z			
				ele.OutPut_Force_2(Forcefilename, step);	//�O�͂̌v�Z			
				ele.Calc_Ohm();								//�d�וۑ����̌v�Z
			}
			else{
				#pragma omp parallel num_threads(Num) if(OMP_calc)
				{
					ele.Calc_Potential_Laplace();
				}
				ele.Calc_E();								//�d�E�v�Z
			}

			//����̏������݁@���Ԃ̏o��
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