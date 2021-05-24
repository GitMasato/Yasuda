#ifndef MAG_CALC_H		//�C���N���[�h�K�[�h
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
	int pos_R;					//�ʒu�ߓ_ 
	int pos_Z;					//�ʒu�ߓ_ 
};

struct OutPut_Coil_info{

	int Coil_Pos_X;				// �R�C���ʒu X ����
	int Coil_Pos_Y;				// �R�C���ʒu Y ����
	int Coil_Pos_Z;				// �R�C���ʒu Z ����
	double coil_pos_x;			// �R�C���ʒu�@X ����
	double coil_pos_y;			// �R�C���ʒu�@Y ����
	double coil_pos_z;			// �R�C���ʒu�@Z ����
	double incline_XZ;			// �R�C���̌X���p�x
	double incline_XY;			// �R�C���̌X���p�x
	double SIN_XY;				// �R�C���̌X���p�x
	double COS_XY;				// �R�C���̌X���p�x
	double SIN_XZ;				// �R�C���̌X���p�x
	double COS_XZ;				// �R�C���̌X���p�x

};

class Magnetic_Calc{

private:

	//��{�p�����[�^
	const double PI;			//�~����
	int wire_turn;				//�R�C��������
	int turn_number_height;		//�R�C����������������
	int nodeNum;				//�ߓ_�̑���	
	int WidthNum;				// r �����ߓ_��
	int HeightNum;				// z �����ߓ_��
	int convergentCount;		//�������Ă��Ȃ��ߓ_��
	int thread;					//���񉻂̃X���b�h��	
	bool OMP_CALC;				//���񉻂��邩
	double ElmSize;				//�v�f����
	double areaWidth;			//�v�Z�̈�R��������
	double areaHeight;			//�v�Z�̈�Z��������
	double permeability;		//�^��̓�����
	double wire_diameter;		//�������a     m
	double wire_totallength;	//��������     m
	double coil_length;			//�R�C������ 
	double coil_diameter;		//�R�C���������~���̒��a�@ m
	double endpoint;			//�v�Z�̈扺�[����R�C����[�܂ł̋����@ m
	double distance_coils;		//�R�C���Ԃ̋���
	double I;					//����d��
	double W;					//����d��    W
	double p;					//�d�C��R���@��m
	double R;					//�R�C���S�̂̒�R ��
	double TFC;					//��������
	double SOR;					//�����W��		

	//�z��f�[�^
	vector<double> Potential;		//�ߓ_�̃x�N�g���|�e���V����
	vector<double> CurrentDensity;	//�ߓ_�̓d�����x
	vector<double> B_r;				//�ߓ_�̎������x
	vector<double> B_z;				//�ߓ_�̎������x
	vector<double> coefficient_1;	//�x�N�g���|�e���V�����v�Z�p�̌W��
	vector<double> coefficient_2;	//�x�N�g���|�e���V�����v�Z�p�̌W��
	vector<double> coefficient_3;	//�x�N�g���|�e���V�����v�Z�p�̌W��
	vector<double> coefficient_4;	//�x�N�g���|�e���V�����v�Z�p�̌W��
	bool * CalcFlag;				//�ߓ_���|�e���V�����v�Z���邩�̃t���O
	vector<Electrode_info> electrode_1;
	vector<Electrode_info> electrode_2;
	vector<OutPut_Coil_info> output_coil;

	//���q�����v�Z�p�̏o�̓p�����[�^
	int nodeNum_output;			//�ߓ_�̑���	
	int WidthNum_output;		// X �����ߓ_��
	int DepthNum_output;		// Y �����ߓ_��
	int HeightNum_output;		// Z �����ߓ_��
	double ElmSize_output;		//�v�f����
	double areaWidth_output;	//�v�Z�̈� X ���������@�@�~�����W�n�́@�R�C�����[���̒��S�ʒu����
	double areaDepth_output;	//�v�Z�̈� Y ���������@
	double areaHeight_output;	//�v�Z�̈� Z ��������	
	vector<double> Output_B_x_0;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_y_0;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_z_0;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_x_1;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_y_1;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_z_1;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_x_2;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_y_2;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_z_2;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_x_3;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_y_3;	//�ߓ_�̎������x�@�~�����W
	vector<double> Output_B_z_3;	//�ߓ_�̎������x�@�~�����W
	

public:

	Magnetic_Calc();							//�R���X�g���N�^
	~Magnetic_Calc();							//�f�X�g���N�^
		
	int Get_thread(){ return thread; }			//�X���b�h����Ԃ�
	bool Get_OMP_Calc(){ return OMP_CALC; }		//���񉻂��邩
	void Calc_CoilPosition();					//�R�C���̈ʒu����
	void Setting_Boundary_Condition();			//���E�����̐ݒ�
	void Calc_Potential_Poisson();				//�|�A�\���������@�x�N�g���|�e���V���������߂�
	void Calc_B();								//�������x�v�Z
	void Cal_B_Cylindrical();				    //�������x�f�[�^���~�����W�n����f�J���g���W�֕ϊ��ɕϊ�
	void Cal_Temp(char *filename_0, char *filename_1, char *filename_2, char *filename_3);

	void OutPut_Coil(char *filename_1, char *filename_2);			//�R�C���ʒu���o�͂���
	void OutPut_Potential(char *filename_1, char *filename_2);		//�x�N�g���|�e���V�����f�[�^���o�͂���
	void OutPut_B_1(char *filename_1, char *filename_2);			//���E�f�[�^���o�͂���
	void OutPut_B_2(char *filename_1, char *filename_2);			//���E�f�[�^���o�͂���
	void OutPut_B(char *filename_1, char *filename_2);				//���E�f�[�^���o�͂���

	void OutPut_B_Cylindrical(char *filename_0_a, char *filename_0_b,
		char *filename_1_a, char *filename_1_b, 
		char *filename_2_a, char *filename_2_b,
		char *filename_3_a, char *filename_3_b);	//�f�J���g���W�n�̎������x���������o�͂���Ă��邩�m�F	
	
	
	void OutPut_MagFieldBinary(char *filename_0, char *filename_1, char *filename_2, char *filename_3);	//�����v�Z�p�̎��E�f�[�^���o�͂���							
	void OutPut_Condition(char *filename);							//�v�Z�������o�͂���

};
#endif