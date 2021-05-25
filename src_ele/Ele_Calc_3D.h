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

struct Ptcl_info{

	vector<int> contain;
	vector<int> node;						
	vector<int> ele;		
	double Ptcl_Pos_X;		//���q�̈ʒu
	double Ptcl_Pos_Y;		//���q�̈ʒu
	double Ptcl_Pos_Z;		//���q�̈ʒu
	double Ptcl_Radius;		//���q�̔��a
	double Ptcl_Charge;		//���q�̑ѓd��
	double sp_G;			//��d	
};

struct Cylinder_info{

	int pos_X;						//�ʒu�ߓ_ 
	int pos_Y;						//�ʒu�ߓ_ 
	int pos_Z;						//�ʒu�ߓ_ 
	int radius;
	int length;
	int diele_radius;
	double Pos_X;					//��[�_�̍��W
	double Pos_Y;					//��[�_�̍��W
	double Pos_Z;					//��[�_�̍��W
	double Radius;					//���a
	double Length;					//����
	double Dielectric_Radius;		//�U�d��
	vector<int> cyl_mesh;	    	//�ʒu�ߓ_ 
};

struct Hollow_Cylinder_info{

	int pos_X;						//�ʒu�ߓ_ 
	int pos_Y;						//�ʒu�ߓ_ 
	int pos_Z;						//�ʒu�ߓ_ 
	int inside_radius;
	int outside_radius;
	int length;
	double Pos_X;					//�[�_�̍��W
	double Pos_Y;					//�[�_�̍��W
	double Pos_Z;					//�[�_�̍��W
	double Inside_Radius;			//���a
	double Outside_Radius;			//���a
	double Length;					//����
	vector<int> hollow_cyl_mesh;	//�ʒu�ߓ_ 
};

struct Wall_info{

	int pos_X[4];					//�ʒu�ߓ_ 
	int pos_Y[4];					//�ʒu�ߓ_ 
	int pos_Z[4];					//�ʒu�ߓ_ 
	int length;						//����
	double Pos_X[4];				//�[�_�̍��W
	double Pos_Y[4];				//�[�_�̍��W
	double Pos_Z[4];				//�[�_�̍��W
	double Length;					//����
	vector<int> wall_mesh;			//�ʒu�ߓ_ 
};

class Electric_Calc{

private:
		
	//��{�p�����[�^
	int Type_boundary_X;						//���E�����̃p�^�[��
	int Type_boundary_Y;						//���E�����̃p�^�[��
	int cellNum_2D;								//�ߓ_�̑���
	int cellNum_3D;								//�ߓ_�̑���	
	int WidthNum;								//x�����ߓ_��
	int DepthNum;								//y�����ߓ_��
	int HeightNum;								//z�����ߓ_��	
	int convergentCount;						//�������Ă��Ȃ��ߓ_��
	int thread;									//���񉻂̃X���b�h��
	bool OMP_CALC;								//���񉻂��邩	
	bool unsteady_calc;							//����̌v�Z��
	bool load_ptcl_Poisson;						//���q�f�[�^��ǂݍ��ނ��@�|�A�\����
	const double PI;							//�~����

	double ElmSize;								//�v�f����
	double areaWidth;							//�v�Z�̈�X��������
	double areaDepth;							//�v�Z�̈�Y��������
	double areaHeight;							//�v�Z�̈�Z��������
	double spElePerm;							//�^��̗U�d��
	double AirRelElePerm;						//��C��U�d��
	double Volt;								//����d��	
	double TFC;									//��������
	double SOR;									//�����W��

	vector<int> bcc;					        //�v�Z����ߓ_	(�d�ʓ���p)
	vector<int> fcc;						    //�v�Z����ߓ_
	vector<int> cbc;					        //�v�Z����ߓ_
	vector<int> cfc;						    //�v�Z����ߓ_
	vector<int> ccb;							//�v�Z����ߓ_
	vector<int> ccf;							//�v�Z����ߓ_
	vector<int> bbb;					        //�v�Z����ߓ_	(�U�d������p)
	vector<int> fbb;						    //�v�Z����ߓ_
	vector<int> bfb;					        //�v�Z����ߓ_
	vector<int> ffb;						    //�v�Z����ߓ_
	vector<int> bbf;					        //�v�Z����ߓ_	
	vector<int> fbf;						    //�v�Z����ߓ_
	vector<int> bff;					        //�v�Z����ߓ_
	vector<int> fff;						    //�v�Z����ߓ_
	vector<double> Potential;					//�ߓ_�̓d��
	vector<double> E_x;							//�ߓ_�̓d�E
	vector<double> E_y;							//�ߓ_�̓d�E
	vector<double> E_z;							//�ߓ_�̓d�E
	vector<double> ElePerm;						//�v�f�̗U�d��
	bool * CalcFlag;							//�ߓ_���|�e���V�����v�Z���邩�̃t���O

	//����d�E�v�Z���s���ꍇ�̂�
	int total_step;								//���X�e�b�v
	int phase;									//������
	int WriteAVS;								//AVS����̕�����
	int phase_No;								//�v���X�̓d��				
	bool phase_change_flag;						//�d�E��؂�ւ��邩
	double span;								//�ꎞ�ϐ�
	double delta_T;								//���U������
	double total_T;								//�v�Z����
	double frequency;							//���g��
	double rise_time;							//�����オ�莞��
		
	//���q�f�[�^	�|�A�\���̎��������ꍇ�̂�	
	int ptclNum;								//���q��
	double conductivity;						//���q�̓��d��
	double ParticleRelElePerm;					//���q��U�d��	
	vector<Ptcl_info> P_C;						//Ptcl_info�ɃA�N�Z�X���邽�߂̃C���X�^���X	
	vector<int> ptcl_num_node;					//�ߓ_�ɔz�u����闱�q�̔ԍ�			
	vector<int> ptcl_num_ele;					//�v�f�ɔz�u����闱�q�̔ԍ�	
	vector<double> Conductivity;				//�v�f�̓��d��	
	vector<double> Charge;						//�ߓ_�̓d��
	
	//�~��������ꍇ
	int cylinder_Xdir_Num;						
	int cylinder_Ydir_Num;						
	vector<Cylinder_info> cylinder_Xdir;
	vector<Cylinder_info> cylinder_Ydir;
	double CylinderElePerm;

	//����~��������ꍇ
	int hollow_cylNum;
	vector<Hollow_Cylinder_info> hollow_cylinder;
	double hollow_cylinderElePerm;

	//�����̂�����ꍇ
	int wallNum;
	vector<Wall_info> wall;
	double wallElePerm;
		
public:

	Electric_Calc();												//�R���X�g���N�^
	~Electric_Calc();												//�f�X�g���N�^

	int Get_AVS(){ return WriteAVS; }								//�X���b�h����Ԃ�
	int Get_thread(){ return thread; }								//�X���b�h����Ԃ�
	bool Get_OMP_Calc(){ return OMP_CALC; }							//���񉻂��邩
	bool Get_unsteady_calc(){ return unsteady_calc; }				//����̌v�Z��
	bool Get_load_ptcl_Poisson(){ return load_ptcl_Poisson; }		//���q���l�����邩�E�|�A�\���̎���
	int Return_step(){ return total_step; }							//�v�Z�X�e�b�v����Ԃ��֐�

	void Setting_Boundary_Condition_V();							//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_X_n();						//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_X_l();						//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_Y_n();						//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_Y_l();						//���E�����̐ݒ�
	void Setting_Boundary_Condition_V_Z_n();						//���E�����̐ݒ�

	void Setting_Boundary_Condition_e();							//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_X_n_Y_n_Z_n();				//���E�����̐ݒ�

	void Setting_Boundary_Condition_e_X_n_corner();					//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_Y_n_corner();					//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_Z_n_corner();					//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_X_n_Y_n_corner();				//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_X_n_Z_n_corner();				//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_Y_n_Z_n_corner();				//���E�����̐ݒ�
	void Setting_Boundary_Condition_e_X_n_Y_n_Z_n_corner();			//���E�����̐ݒ�
	
	void Setting_Boundary_Condition_1();							//���E�����̐ݒ�
	void Setting_Boundary_Condition_2(int step, char *filename);	//���E�����̐ݒ�
	void Ptcl_Data(char *filename);									//���q�f�[�^�̓ǂݍ���		
	void Particle_Node_Mapping();									//���q��ߓ_��ɐݒu
	void Particle_Ele_Mapping();									//���q��v�f��ɐݒu
	void Particle_count_charge();									//���q���ƂɁA�d�ׂ�z�u����ߓ_���𐔂���
	void Particle_charge_Mapping();									//�d�ׂ�z�u����
	void Particle_clear();											//�s�v�ȃ��������������
	void Calc_Potential_Laplace();									//���v���X������		
	void Calc_Potential_Poisson();									//�|�A�\��������
	void Calc_E();													//�d�E�v�Z
	void Calc_E_X_n();												//�d�E�v�Z
	void Calc_E_X_l();												//�d�E�v�Z
	void Calc_E_Y_n();												//�d�E�v�Z
	void Calc_E_Y_l();												//�d�E�v�Z
	void Calc_E_Z_n();												//�d�E�v�Z

	void Calc_Ohm();												//�I�[���̎��v�Z		
	void Ini_P();													//�d�ʕ��z��������
	void OutPut_Potential(char *filename_1, char *filename_2);		//�d�ʃf�[�^���o�͂���
	void OutPut_E(char *filename_1, char *filename_2);				//�d�E�f�[�^���o�͂���
	void OutPut_Charge(char *filename_1, char *filename_2);			//�d�׃f�[�^���o�͂���
	void OutPut_Ptcl(char *filename);								//���q�ʒu���o�͂���
	void OutPut_Force_1(char *filename);							//���q�ɉ����O�͂��o�͂���
	void OutPut_Force_2(char *filename, int step);					//���q�ɉ����O�͂��o�͂���
	void OutPut_EleFieldBinary(char *filename);						//�����v�Z�p�̓d�E�f�[�^���o�͂���
	void GenerateMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat);	//����f�[�^���쐬		
	void WriteMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat, int T_step, int Count);	//����f�[�^����������
	void OutPut_Condition(char *filename);

	int int_string(const string& str){			//string�^��int�^�֕ϊ�����֐�

		int t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	double double_string(const string& str){			//string�^��double�^�֕ϊ�����֐�

		double t;
		stringstream ss;
		ss << str;
		ss >> t;
		return t;
	}

	list<string> split(string str, string delim){

		list<string> result;
		int cutAt;

		while ((cutAt = int(str.find_first_of(delim))) != str.npos){

			if (cutAt > 0){

				result.push_back(str.substr(0, cutAt));
			}

			str = str.substr(cutAt + 1);

		}

		if (str.length() > 0){

			result.push_back(str);
		}

		return result;
	}
};
#endif