#ifndef ELE_CALC_H		//�C���N���[�h�K�[�h
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

	int pos_X;						//�ʒu�ߓ_ 
	int pos_Y;						//�ʒu�ߓ_ 
	int pos_Z;						//�ʒu�ߓ_ 
	int radius;
	int diele_radius;
	double Pos_X;					//��[�_�̍��W
	double Pos_Y;					//��[�_�̍��W
	double Pos_Z;					//��[�_�̍��W
	double Radius;					//���a
	double Dielectric_Radius;		//�U�d��
};

struct Wall_info{

	int pos_X[4];					//�ʒu�ߓ_ 
	int pos_Y[4];					//�ʒu�ߓ_ 
	int pos_Z[4];					//�ʒu�ߓ_ 
	double Pos_X[4];				//�[�_�̍��W
	double Pos_Y[4];				//�[�_�̍��W
	double Pos_Z[4];				//�[�_�̍��W
};

class Electric_Calc{

private:

	//��{�p�����[�^
	const double PI;							//�~����
	int Type_boundary_X;						//���E�����̃p�^�[��
	int cellNum;								//�ߓ_�̑���	
	int WidthNum;								//x�����ߓ_��
	int HeightNum;								//z�����ߓ_��
	int convergentCount;						//�������Ă��Ȃ��ߓ_��
	int thread;									//���񉻂̃X���b�h��	
	bool OMP_CALC;								//���񉻂��邩
	double ElmSize;								//�v�f����
	double areaWidth;							//�v�Z�̈�X��������
	double areaHeight;							//�v�Z�̈�Z��������
	double spElePerm;							//�^��̗U�d���D
	double AirRelElePerm;						//��C��U�d��
	double Volt;								//����d��
	double TFC;									//��������
	double SOR;									//�����W��

	vector<int> bcc;					        //�v�Z����ߓ_�̍��ߓ_(�d�ʓ���p)
	vector<int> fcc;						    //�v�Z����ߓ_�̉E�ߓ_
	vector<int> ccb;							//�v�Z����ߓ_�̏�ߓ_
	vector<int> ccf;							//�v�Z����ߓ_�̉��ߓ_
	vector<int> bcb;							//�v�Z����ߓ_�̍����ߓ_(�U�d������p)
	vector<int> bcf;							//�v�Z����ߓ_�̉E���ߓ_
	vector<int> fcb;							//�v�Z����ߓ_�̍���ߓ_
	vector<int> fcf;							//�v�Z����ߓ_�̉E��ߓ_

	vector<double> ElePerm;						//�v�f�̗U�d��
	vector<double> Potential;					//�ߓ_�̓d��
	vector<double> E_x;							//�ߓ_�̓d�E
	vector<double> E_z;							//�ߓ_�̓d�E
	bool * CalcFlag;							//�ߓ_���|�e���V�����v�Z���邩�̃t���O

	//�~������ꍇ
	int cylinder_Num;
	vector<Cylinder_info> cylinder;
	double CylinderElePerm;
	
	//�l�p������ꍇ
	int wallNum_electrode;
	vector<Wall_info> wall_electrode;
	int wallNum_dielectric;
	vector<Wall_info> wall_dielectric;
	double wallElePerm;

public:

	Electric_Calc();												//�R���X�g���N�^
	~Electric_Calc();												//�f�X�g���N�^

	int Get_thread(){ return thread; }								//�X���b�h����Ԃ�
	bool Get_OMP_Calc(){ return OMP_CALC; }							//���񉻂��邩
	void Setting_Boundary_Condition_V();							//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_X_n();						//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_X_r();						//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_Z_n();						//���E�����̐ݒ�

	void Setting_Boundary_Condition_e();							//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_X_n_Z_n();					//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_X_r_Z_n();					//���E�����̐ݒ�

	void Setting_Boundary_Condition();								//���E�����̐ݒ�
	void Calc_Potential_Laplace();									//���v���X������
	void Calc_E();													//�d�E�v�Z
	void Calc_E_X_n();												//�d�E�v�Z
	void Calc_E_X_r();												//�d�E�v�Z
	void Calc_E_Z_n();												//�d�E�v�Z

	void OutPut_Potential(char *filename_1, char *filename_2);		//�d�ʃf�[�^���o�͂���
	void OutPut_E(char *filename_1, char *filename_2);				//�d�E�f�[�^���o�͂���
	void OutPut_EleFieldBinary(char *filename);						//�����v�Z�p�̓d�E�f�[�^���o�͂���							
	void OutPut_Condition(char *filename);
	
};
#endif