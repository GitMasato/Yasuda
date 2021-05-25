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

//main�֐���
int main(){
	
	std::ios::sync_with_stdio(false);
	
	//exe���̎��s�t�@�C���Ɠ����t�H���_��OUTPUT�t�H���_���쐬���Ă���
	//�o�̓t�@�C���̖��O
	char E_Field_1[64] = "1FieldData_2D.bin";						//���v���X��̃f�[�^
	char End_Laplace_Potentialfilename_1[64] = "P_Laplace.fld";		//���v���X��̓d�ʃf�[�^
	char End_Laplace_Potentialfilename_2[64] = "P_Laplace.dat";		//���v���X��̓d�ʃf�[�^
	char End_Laplace_E_filename_1[64] = "E_Laplace.fld";			//�|�A�\����̓d�ʃf�[�^
	char End_Laplace_E_filename_2[64] = "E_Laplace.dat";			//�|�A�\����̓d�ʃf�[�^
	char memo[64] = "memo.txt";
		
	//�ϐ��̏�������
	Electric_Calc ele;
	Profiler p;
	int Num = ele.Get_thread();
	bool OMP_calc = ele.Get_OMP_Calc();

	//���E�����̐ݒ�
	ele.Setting_Boundary_Condition_V();
	ele.Setting_Boundary_Condition_e();
	ele.Setting_Boundary_Condition();

	//���v���X�̎�
	#pragma omp parallel num_threads(Num) if(OMP_calc)
	{
		ele.Calc_Potential_Laplace();
	}

	//���v���X�̎�����������̓d�ʕ��z���o��
	ele.OutPut_Potential(End_Laplace_Potentialfilename_1, End_Laplace_Potentialfilename_2);

	//���v���X�v�Z��̓d�E�v�Z
	ele.Calc_E();

	//���v���X�̎�����������̓d�E���z���o��
	ele.OutPut_E(End_Laplace_E_filename_1, End_Laplace_E_filename_2);
	
	//�����v�Z�p�̓d�E�f�[�^���o��
	ele.OutPut_EleFieldBinary(E_Field_1);

	//�v�Z���������o�͂��Ă���
	ele.OutPut_Condition(memo);

	//�v�Z���Ԃ��o�͂���
	ofstream fout;
	fout.open("2_elapsed_time.txt", ios_base::out | ios_base::trunc);
	fout << "erapsed time = " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
	fout.close();

}