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

//main�֐���
int main(){

	std::ios::sync_with_stdio(false);
	Profiler p;

	//exe���̎��s�t�@�C���Ɠ����t�H���_��OUTPUT�t�H���_���쐬���Ă���
	//�o�̓t�@�C���̖��O
	char B_Field_0[64] = "1Mag_FieldData_3D_0.bin";				//�|�A�\����̗��q�����v�Z�p�f�[�^
	char B_Field_1[64] = "1Mag_FieldData_3D_1.bin";				//�|�A�\����̗��q�����v�Z�p�f�[�^
	char B_Field_2[64] = "1Mag_FieldData_3D_2.bin";				//�|�A�\����̗��q�����v�Z�p�f�[�^
	char B_Field_3[64] = "1Mag_FieldData_3D_3.bin";				//�|�A�\����̗��q�����v�Z�p�f�[�^

	char Cartesian_0_a[64] = "Cartesian_0.fld";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^
	char Cartesian_0_b[64] = "Cartesian_0.dat";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^
	char Cartesian_1_a[64] = "Cartesian_1.fld";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^
	char Cartesian_1_b[64] = "Cartesian_1.dat";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^
	char Cartesian_2_a[64] = "Cartesian_2.fld";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^
	char Cartesian_2_b[64] = "Cartesian_2.dat";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^
	char Cartesian_3_a[64] = "Cartesian_3.fld";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^
	char Cartesian_3_b[64] = "Cartesian_3.dat";						//�|�A�\����̗��q�����v�Z�p�f�[�^���`�F�b�N���邽�߂̃f�[�^

	char Poisson_Potentialfilename_1[64] = "P_Poisson.fld";		//�|�A�\����̃x�N�g���|�e���V�����f�[�^
	char Poisson_Potentialfilename_2[64] = "P_Poisson.dat";		//�|�A�\����̃x�N�g���|�e���V�����f�[�^
	char Poisson_B_filename_1[64] = "B_Poisson_1.fld";			//�|�A�\����̎������x�f�[�^
	char Poisson_B_filename_2[64] = "B_Poisson_1.dat";			//�|�A�\����̎������x�f�[�^
	char Poisson_B_filename_3[64] = "B_Poisson_2.fld";			//�|�A�\����̎������x�f�[�^
	char Poisson_B_filename_4[64] = "B_Poisson_2.dat";			//�|�A�\����̎������x�f�[�^
	char Poisson_Br[64] = "Br.csv";                             //�R�C���̔����̈ʒu�̎������x�f�[�^
	char Poisson_Bz[64] = "Bz.csv";								//���S����̎������x�f�[�^
	
	char Coil_filename_1[64] = "Coil_Position.fld";				//�R�C���ʒu�f�[�^
	char Coil_filename_2[64] = "Coil_Position.dat";				//�R�C���f�[�^
	
	char Peidong_Bz_0[64] = "Bz_dis_0.csv";						//�ꎞ�I�ȃf�[�^
	char Peidong_Bz_1[64] = "Bz_dis_1.csv";						//�ꎞ�I�ȃf�[�^
	char Peidong_Bz_2[64] = "Bz_dis_2.csv";						//�ꎞ�I�ȃf�[�^
	char Peidong_Bz_3[64] = "Bz_dis_3.csv";						//�ꎞ�I�ȃf�[�^
	char memo[64] = "memo.txt";
	
	
	
	//�ϐ��̏�������
	Magnetic_Calc mag;
	int Num = mag.Get_thread();
	bool OMP_calc = mag.Get_OMP_Calc();

	//�R�C��1���̈ʒu�Ȃǌv�Z
	mag.Calc_CoilPosition();

	//���E�����̐ݒ�
	#pragma omp parallel num_threads(Num) if(OMP_calc)
	{
		mag.Setting_Boundary_Condition();
	}

	//�R�C���ʒu���o��
	mag.OutPut_Coil(Coil_filename_1, Coil_filename_2);	
	
	//�|�A�\���̎�
	#pragma omp parallel num_threads(Num) if(OMP_calc)
	{
		mag.Calc_Potential_Poisson();
	}

	//�|�A�\���̎�����������̃x�N�g���|�e���V�������z���o��
	mag.OutPut_Potential(Poisson_Potentialfilename_1, Poisson_Potentialfilename_2);
	
	//�������x�v�Z
	mag.Calc_B();				

	//�|�A�\���̎�����������̎������x���z���o��
	mag.OutPut_B_1(Poisson_B_filename_1, Poisson_B_filename_2);
	mag.OutPut_B_2(Poisson_B_filename_3, Poisson_B_filename_4);
	mag.OutPut_B(Poisson_Br,Poisson_Bz);

	//�������x�f�[�^���~�����W�n����f�J���g���W�֕ϊ��ɕϊ�
	mag.Cal_B_Cylindrical();
	mag.Cal_Temp(Peidong_Bz_0, Peidong_Bz_1, Peidong_Bz_2, Peidong_Bz_3);

	//�f�J���g���W�n�̎������x���������o�͂���Ă��邩�m�F
	mag.OutPut_B_Cylindrical(Cartesian_0_a, Cartesian_0_b, 
		Cartesian_1_a, Cartesian_1_b, 
		Cartesian_2_a, Cartesian_2_b, 
		Cartesian_3_a, Cartesian_3_b);

	//�����v�Z�p�̓d�E�f�[�^���o��
	mag.OutPut_MagFieldBinary(B_Field_0, B_Field_1, B_Field_2, B_Field_3);

	//�v�Z���������o��
	mag.OutPut_Condition(memo);
	
	//�v�Z���Ԃ��o�͂���
	ofstream fout;
	fout.open("2_elapsed_time.txt", ios_base::out | ios_base::trunc);
	fout << "erapsed time = " << p.print() / 3600 << ":" << (p.print() % 3600) / 60 << ":" << (p.print() % 3600) % 60 << endl;
	fout.close();
	
}