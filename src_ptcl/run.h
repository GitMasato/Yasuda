#ifndef RUN_H		//�C���N���[�h�K�[�h
#define RUN_H		//�C���N���[�h�K�[�h

#include <iostream>
#include <iomanip>
#include <numeric>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <new>
#include <random>
#include <algorithm>
#include "DataIO.h"

using namespace std;

class Run
{

private:

	DataIO d;
	bool load_ptcl_number_diameter;
	int ptclNum;
	double min_Pos[3], max_Pos[3];
	double min_Velo[3], max_Velo[3];
	double min_omega[3], max_omega[3];
	double diameter_average, diameter_stadard_deviation;
	double ave_sp_gravity, standev_sp_gravity;
	double ave_Charge, standev_Charge;
	double ave_ele_conductivity, standev_ele_conductivity;
	double ave_permittivity, standev_permittivity;
	double ave_permeability, standev_permeability;
	double ave_adhesion, standev_adhesion;
	double ave_rolling_fri, standev_rolling_fri;
	double PI;
	vector<double> Pos_X, Pos_Y, Pos_Z;
	vector<double> Velo_X, Velo_Y, Velo_Z;
	vector<double> Omega_X, Omega_Y, Omega_Z;
	vector<double> Diameter;
	vector<double> SP_Gravity;
	vector<double> Charge;
	vector<double> Conductivity;
	vector<double> Permittivity;
	vector<double> Permeability;
	vector<double> Adhesion;
	vector<double> Rolling_Fri;

public:

	Run();
	~Run();

	void Load(char *filename_1, char *filename_2);			//  �ݒ�t�@�C����ǂݍ���
	void Particle_Initial_Placement();						//	���q���ɒl��ݒ�

	void Sp_gravity();										//�@��d��ݒ�
	void SizeCal();											//  ���a�ݒ�
	void ChargeCal();										//  ��ѓd�ʐݒ�
	void EleconductCal();									//	���d���ݒ�
	void PermitCal();										//	�U�d���ݒ�
	void PermeaCal();										//	�������ݒ�
	void AdhesionCal();										//	�t���͌W���ݒ�
	void RollingFriction();									//	�]���薀�C�W���ݒ�
	void PlaceCal();										//  ���q�ʒu�ݒ�
	void VeloCal();											//  ���q���x�ݒ�
	void OmegaCal();										//  ���q�p���x�ݒ�
	void Output_ptcl(char *filename);						//	�f�[�^�o��

};

#endif
//�C���N���[�h�K�[�h

