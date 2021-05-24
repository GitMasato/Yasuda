#include "Mag_Calc.h"

Magnetic_Calc::Magnetic_Calc() :PI(3.14159265358979323846264338328){		//�R���X�g���N�^�D�ϐ��̏����������Ă���D

	cout << "�ϐ��E�v�Z�����̐ݒ�" << endl;

	//�ϐ��̐ݒ�
	ElmSize = 100.0e-6;				//�v�Z���b�V�� m
	areaWidth = 200.0e-3;			//�v�Z�͈� m�@r����
	areaHeight = 616.3e-3;			//�v�Z�͈� m�@z����

	endpoint = 200.0e-3;			//�v�Z�̈扺�[����R�C����[�܂ł̋����@ m
	distance_coils = 140.3e-3;		//�R�C���\�ʊԂ̋����@ m
	wire_diameter = 500.0e-6;		//�������a m
	coil_length = 38.0e-3;			//�R�C������ m
	coil_diameter = 10.0e-3;		//�R�C���������~���̒��a m
	wire_turn = 5000;				//����������	
	I = 2.0;						//�d��        A
	p = 1.724e-8;					//�d�C��R���@��m
	permeability = 1.25633e-6;		//�^��̓�����  H/m
	OMP_CALC = true;				//���񉻂��邩
	thread = 4;						//����
	TFC = 1.0e-9;					//��������
	SOR = 1.8;						//�����W��	
	
	//�����v�Z�p�̃f�[�^�o�͗p�̕ϐ��ݒ�
	output_coil.clear();
	output_coil.resize(4);

	ElmSize_output = 500.0e-6;		//�v�f����
	areaWidth_output = 50.0e-3;		//�v�Z�̈� X ��������
	areaDepth_output = 50.0e-3;		//�v�Z�̈� Y ���������@
	areaHeight_output = 50.0e-3;	//�v�Z�̈� Z ��������

	output_coil[0].incline_XY = 45 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[0].incline_XZ = 35.26 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[0].coil_pos_x = -37.5e-3;			//�R�C���ʒu�@X ����   ��[�̒��S����
	output_coil[0].coil_pos_y = -37.5e-3;			//�R�C���ʒu�@Y ����   ��[�̒��S����
	output_coil[0].coil_pos_z = -37.5e-3;			//�R�C���ʒu�@Z ����   ��[�̒��S����	

	output_coil[1].incline_XY = 135 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[1].incline_XZ = 35.26 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[1].coil_pos_x = 87.5e-3;			//�R�C���ʒu�@X ����   ��[�̒��S����
	output_coil[1].coil_pos_y = -37.5e-3;			//�R�C���ʒu�@Y ����   ��[�̒��S����
	output_coil[1].coil_pos_z = -37.5e-3;			//�R�C���ʒu�@Z ����   ��[�̒��S����	

	output_coil[2].incline_XY = 225 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[2].incline_XZ = 35.26 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[2].coil_pos_x = 87.5e-3;			//�R�C���ʒu�@X ����   ��[�̒��S����
	output_coil[2].coil_pos_y = 87.5e-3;			//�R�C���ʒu�@Y ����   ��[�̒��S����
	output_coil[2].coil_pos_z = -37.5e-3;			//�R�C���ʒu�@Z ����   ��[�̒��S����	

	output_coil[3].incline_XY = 315 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[3].incline_XZ = 35.26 * (PI / 180);	// �R�C���̌X���p�x (���W�A���ɕϊ�)
	output_coil[3].coil_pos_x = -37.5e-3;			//�R�C���ʒu�@X ����   ��[�̒��S����
	output_coil[3].coil_pos_y = 87.5e-3;			//�R�C���ʒu�@Y ����   ��[�̒��S����
	output_coil[3].coil_pos_z = -37.5e-3;			//�R�C���ʒu�@Z ����   ��[�̒��S����	


	//�ȍ~�͕ϐ��̏������E�p�����[�^�̗��U��
	output_coil[0].Coil_Pos_X = (int)(output_coil[0].coil_pos_x / ElmSize_output + 0.5);
	output_coil[0].Coil_Pos_Y = (int)(output_coil[0].coil_pos_y / ElmSize_output + 0.5);
	output_coil[0].Coil_Pos_Z = (int)(output_coil[0].coil_pos_z / ElmSize_output + 0.5);
	output_coil[1].Coil_Pos_X = (int)(output_coil[1].coil_pos_x / ElmSize_output + 0.5);
	output_coil[1].Coil_Pos_Y = (int)(output_coil[1].coil_pos_y / ElmSize_output + 0.5);
	output_coil[1].Coil_Pos_Z = (int)(output_coil[1].coil_pos_z / ElmSize_output + 0.5);
	output_coil[2].Coil_Pos_X = (int)(output_coil[2].coil_pos_x / ElmSize_output + 0.5);
	output_coil[2].Coil_Pos_Y = (int)(output_coil[2].coil_pos_y / ElmSize_output + 0.5);
	output_coil[2].Coil_Pos_Z = (int)(output_coil[2].coil_pos_z / ElmSize_output + 0.5);
	output_coil[3].Coil_Pos_X = (int)(output_coil[3].coil_pos_x / ElmSize_output + 0.5);
	output_coil[3].Coil_Pos_Y = (int)(output_coil[3].coil_pos_y / ElmSize_output + 0.5);
	output_coil[3].Coil_Pos_Z = (int)(output_coil[3].coil_pos_z / ElmSize_output + 0.5);

	output_coil[0].SIN_XY = sin(output_coil[0].incline_XY);
	output_coil[0].COS_XY = cos(output_coil[0].incline_XY);
	output_coil[1].SIN_XY = sin(output_coil[1].incline_XY);
	output_coil[1].COS_XY = cos(output_coil[1].incline_XY);
	output_coil[2].SIN_XY = sin(output_coil[2].incline_XY);
	output_coil[2].COS_XY = cos(output_coil[2].incline_XY);
	output_coil[3].SIN_XY = sin(output_coil[3].incline_XY);
	output_coil[3].COS_XY = cos(output_coil[3].incline_XY);

	output_coil[0].SIN_XZ = sin(output_coil[0].incline_XZ);
	output_coil[0].COS_XZ = cos(output_coil[0].incline_XZ);
	output_coil[1].SIN_XZ = sin(output_coil[1].incline_XZ);
	output_coil[1].COS_XZ = cos(output_coil[1].incline_XZ);
	output_coil[2].SIN_XZ = sin(output_coil[2].incline_XZ);
	output_coil[2].COS_XZ = cos(output_coil[2].incline_XZ);
	output_coil[3].SIN_XZ = sin(output_coil[3].incline_XZ);
	output_coil[3].COS_XZ = cos(output_coil[3].incline_XZ);	

	WidthNum_output = (int)(areaWidth_output / ElmSize_output + 0.5) +1;
	DepthNum_output = (int)(areaDepth_output / ElmSize_output + 0.5) +1;
	HeightNum_output = (int)(areaHeight_output / ElmSize_output + 0.5) +1;
	nodeNum_output = WidthNum_output * DepthNum_output * HeightNum_output;
	if (nodeNum_output == 0)	nodeNum_output = 1;
	
	Output_B_x_0.assign(nodeNum_output, 0.0);
	Output_B_y_0.assign(nodeNum_output, 0.0);
	Output_B_z_0.assign(nodeNum_output, 0.0);
	Output_B_x_1.assign(nodeNum_output, 0.0);
	Output_B_y_1.assign(nodeNum_output, 0.0);
	Output_B_z_1.assign(nodeNum_output, 0.0);
	Output_B_x_2.assign(nodeNum_output, 0.0);
	Output_B_y_2.assign(nodeNum_output, 0.0);
	Output_B_z_2.assign(nodeNum_output, 0.0);
	Output_B_x_3.assign(nodeNum_output, 0.0);
	Output_B_y_3.assign(nodeNum_output, 0.0);
	Output_B_z_3.assign(nodeNum_output, 0.0);
	
	WidthNum = (int)(areaWidth / ElmSize + 0.5);
	HeightNum = (int)(areaHeight / ElmSize + 0.5);
	if (WidthNum == 0)	WidthNum = 1;

	coefficient_1.assign(WidthNum, 0.0);
	coefficient_2.assign(WidthNum, 0.0);
	coefficient_3.assign(WidthNum, 0.0);
	coefficient_4.assign(WidthNum, 0.0);
	
	for (int i = 0; i < WidthNum; i++){
		coefficient_1[i] = 1.0 / ((4.0 * double(i * ElmSize) * double(i * ElmSize)) + (ElmSize * ElmSize));
		coefficient_2[i] = double(i * ElmSize) * double(i * ElmSize);
		coefficient_3[i] = double(i * ElmSize) * ElmSize * 0.5;
		coefficient_4[i] = double(i * ElmSize) * double(i * ElmSize) * ElmSize * ElmSize;
	}
	
	nodeNum = WidthNum * HeightNum;
	if (nodeNum == 0)	nodeNum = 1;

	Potential.assign(nodeNum, 0.0);
	CurrentDensity.assign(nodeNum, 0.0);
	B_r.assign(nodeNum, 0.0);
	B_z.assign(nodeNum, 0.0);
	CalcFlag = new bool[nodeNum];
	for (int i = 0; i < nodeNum; i++)	CalcFlag[i] = true;

	turn_number_height = (int)(coil_length / wire_diameter);	//�R�C����������������
	wire_totallength = 0.0;										//��������    m		 (��Ōv�Z����̂ŁC�Ƃ肠����0)
	W = 0.0;													//����d��    W		 (��Ōv�Z����̂ŁC�Ƃ肠����0)
	R = 0.0;													//�R�C���S�̂̒�R ��	(��Ōv�Z����̂ŁC�Ƃ肠����0)
	convergentCount = 0;										//�������Ă��Ȃ��ߓ_��	(��Ōv�Z����̂ŁC�Ƃ肠����0)
	electrode_1.clear();
	electrode_1.resize(wire_turn);								//�R�C��1�������̍\���̂����������쐬
	electrode_2.clear();
	electrode_2.resize(wire_turn);								//�R�C��1�������̍\���̂����������쐬


}

Magnetic_Calc::~Magnetic_Calc(){

}

void Magnetic_Calc::Calc_CoilPosition(){		//�R�C���̈ʒu����

	cout << "�R�C��1���̈ʒu���v�Z" << endl;

	//�e�R�C���̈ʒu��ݒ�
	int count_1 = 0, count_2 = 0;
	wire_totallength = 0.0;						//�R�C�������S�̂̒���	

	for (int i = 0; i < wire_turn; i++){

		if (count_1 >= turn_number_height){
			count_1 = 0;
			count_2++;
		}
				
		const double coil_r = (coil_diameter * 0.5) + (wire_diameter * 0.5) + (count_2 * wire_diameter);
		const double coil_z = endpoint + (wire_diameter * 0.5) + (count_1 * wire_diameter);
		
		wire_totallength += coil_r * 2.0 * PI;
		electrode_1[i].pos_R = (int)(coil_r / ElmSize + 0.5);
		electrode_1[i].pos_Z = (int)(coil_z / ElmSize + 0.5);

		electrode_2[i].pos_R = (int)(coil_r / ElmSize + 0.5);
		electrode_2[i].pos_Z = (int)((coil_z + coil_length + distance_coils) / ElmSize + 0.5);

		count_1++;		
	}

	//����d�͌v�Z
	R = p*wire_totallength / (PI * pow(wire_diameter*0.5, 2.0));
	W = I*I*R;

}

void Magnetic_Calc::Setting_Boundary_Condition(){		//���E�����ݒ�

	#pragma omp single
	cout << "�R�C���ɗ����d���̋��E������ݒ�" << endl;
	
	//�d�ɂ̓d�����x����ݒ�
	const double current = permeability * I / ((wire_diameter*0.5) * (wire_diameter*0.5) * PI);
	
	#pragma omp for 
	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum + i;
				
			for (int e = 0; e < electrode_1.size(); e++){

				const double dr = abs(i - electrode_1[e].pos_R) * ElmSize;
				const double dz = abs(k - electrode_1[e].pos_Z) * ElmSize;
				const double distance_temp = (dr * dr) + (dz * dz);
				const double distance = (distance_temp <= 0) ? 0 : sqrt(distance_temp);	
				if (distance <= wire_diameter * 0.5)	CurrentDensity[m] = current;
			}	

			for (int e = 0; e < electrode_2.size(); e++){

				const double dr = abs(i - electrode_2[e].pos_R) * ElmSize;
				const double dz = abs(k - electrode_2[e].pos_Z) * ElmSize;
				const double distance_temp = (dr * dr) + (dz * dz);
				const double distance = (distance_temp <= 0) ? 0 : sqrt(distance_temp);
				if (distance <= wire_diameter * 0.5)	CurrentDensity[m] = current;
			}

		}
	}

	#pragma omp single
	cout << "�[�̋��E�����ݒ�" << endl;

	//���̓d�ʂɊւ��鋫�E�����ݒ�
	#pragma omp for
	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			if (i == 0)					CalcFlag[m] = false;
			if (i == WidthNum-1)		CalcFlag[m] = false;
			if (k == 0)					CalcFlag[m] = false;
			if (k == HeightNum - 1)		CalcFlag[m] = false;
		}
	}
}


/*
void Magnetic_Calc::Setting_Boundary_Condition(){		//���E�����ݒ�

	#pragma omp single
	cout << "�R�C���ɗ����d���̋��E������ݒ�" << endl;

	//�d�ɂ̓d�����x����ݒ�
	const int coil_bottom = (int)(endpoint / ElmSize + 0.5);
	const int coil_upper = coil_bottom + (int)(coil_length / ElmSize + 0.5);
	const int coil_left = (int)((coil_diameter*0.5) / ElmSize + 0.5);
	const int coil_right = coil_left + (int)(wire_diameter * 10.0 / ElmSize + 0.5);
	const double current = permeability * I * 1000.0 / ((coil_upper - coil_bottom) * ElmSize * (coil_right - coil_left) * ElmSize);
	
	#pragma omp for
	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum + i;
			if (((i >= coil_left) && (i <= coil_right)) && ((k >= coil_bottom) && (k <= coil_upper))) 	CurrentDensity[m] = current;
		}
	}

	#pragma omp single
	cout << "�[�̋��E�����ݒ�" << endl;

	//���̓d�ʂɊւ��鋫�E�����ݒ�
	#pragma omp for
	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			if (i == 0)					CalcFlag[m] = false;
			if (i == WidthNum - 1)		CalcFlag[m] = false;
			if (k == 0)					CalcFlag[m] = false;
			if (k == HeightNum - 1)		CalcFlag[m] = false;
		}
	}	
}
*/

void Magnetic_Calc::Calc_Potential_Poisson(){

	int localCount = 0, m;
	double P, Source, P_fcc, P_bcc, P_ccf, P_ccb;

	#pragma omp single
	cout << "�|�A�\���̎����v�Z���܂�" << endl;

	do{
		#pragma omp barrier	
		localCount = 0;

		#pragma omp single	
		convergentCount = 0;

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k = k + 2){
			for (int i = 0; i<WidthNum; i = i + 2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					Source = CurrentDensity[m];
					P_bcc = Potential[m - 1];
					P_fcc = Potential[m + 1];
					P_ccb = Potential[m - WidthNum];
					P_ccf = Potential[m + WidthNum];
					Potential[m] = (1.0 - SOR) * P + SOR * (coefficient_1[i] * ((coefficient_2[i] * (P_fcc + P_bcc + P_ccf + P_ccb)) + (coefficient_3[i] * (P_fcc - P_bcc)) + (coefficient_4[i] * Source)));					

					if (fabs(Potential[m] - P)>TFC)  localCount++;
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k = k + 2){
			for (int i = 1; i<WidthNum; i = i + 2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					Source = CurrentDensity[m];
					P_bcc = Potential[m - 1];
					P_fcc = Potential[m + 1];
					P_ccb = Potential[m - WidthNum];
					P_ccf = Potential[m + WidthNum];
					Potential[m] = (1.0 - SOR) * P + SOR * (coefficient_1[i] * ((coefficient_2[i] * (P_fcc + P_bcc + P_ccf + P_ccb)) + (coefficient_3[i] * (P_fcc - P_bcc)) + (coefficient_4[i] * Source)));

					if (fabs(Potential[m] - P)>TFC)  localCount++;
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k = k + 2){
			for (int i = 1; i<WidthNum; i = i + 2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					Source = CurrentDensity[m];
					P_bcc = Potential[m - 1];
					P_fcc = Potential[m + 1];
					P_ccb = Potential[m - WidthNum];
					P_ccf = Potential[m + WidthNum];
					Potential[m] = (1.0 - SOR) * P + SOR * (coefficient_1[i] * ((coefficient_2[i] * (P_fcc + P_bcc + P_ccf + P_ccb)) + (coefficient_3[i] * (P_fcc - P_bcc)) + (coefficient_4[i] * Source)));

					if (fabs(Potential[m] - P)>TFC)  localCount++;
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k = k + 2){
			for (int i = 0; i<WidthNum; i = i + 2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					Source = CurrentDensity[m];
					P_bcc = Potential[m - 1];
					P_fcc = Potential[m + 1];
					P_ccb = Potential[m - WidthNum];
					P_ccf = Potential[m + WidthNum];
					Potential[m] = (1.0 - SOR) * P + SOR * (coefficient_1[i] * ((coefficient_2[i] * (P_fcc + P_bcc + P_ccf + P_ccb)) + (coefficient_3[i] * (P_fcc - P_bcc)) + (coefficient_4[i] * Source)));

					if (fabs(Potential[m] - P)>TFC)  localCount++;
				}
			}
		}

		#pragma omp critical
		convergentCount += localCount;

		#pragma omp barrier

		#pragma omp single
		cout << convergentCount << "\r";

	} while (convergentCount);
}

void Magnetic_Calc::Calc_B(){

	//�q�����̎������x�v�Z
	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			
			if ((i == 0) || (i == WidthNum - 1) || (k == 0) || (k == HeightNum - 1)){
				B_r[m] = 0;
			}
			else{
				B_r[m] = -(Potential[m + WidthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
			}				
		}
	}

	//�y�����̎������x�v�Z
	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;

			if((i == 0) || (i == WidthNum - 1) || (k == 0) || (k == HeightNum - 1)){
				B_z[m] = 0;
			}
			else {				
				B_z[m] = (Potential[m] / (i*ElmSize)) + (Potential[m + 1] - Potential[m - 1]) / (ElmSize * 2.0);
			}
		}
	}

	//�y�����̎������x�v�Z (���ٓ_�̒��S����)
	for (int k = 0; k < HeightNum; k++){
		
		int m = k * WidthNum;
		B_r[m] = B_r[m + 1];
		B_z[m] = B_z[m + 1];
	}
}

void Magnetic_Calc::Cal_B_Cylindrical(){

	cout << "�������x�f�[�^���~�����W�n����f�J���g���W�֕ϊ�" << endl;
	
	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				const int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				const double temp_x_1 = i * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_y_1 = j * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_z_1 = k * ElmSize_output;		//���W�ϊ���̍��W

				const double temp_x_2 = temp_x_1 * output_coil[0].COS_XY + temp_y_1 * output_coil[0].SIN_XY;	//���W�ϊ���̍��W
				const double temp_y_2 = -temp_x_1 * output_coil[0].SIN_XY + temp_y_1 * output_coil[0].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_2 = temp_z_1;																//���W�ϊ���̍��W
				const double temp_x_3 = temp_x_2 * output_coil[0].COS_XZ + temp_z_2 * output_coil[0].SIN_XZ;	//���W�ϊ���̍��W
				const double temp_y_3 = temp_y_2;																//���W�ϊ���̍��W
				const double temp_z_3 = -temp_x_2 * output_coil[0].SIN_XZ + temp_z_2 * output_coil[0].COS_XZ;	//���W�ϊ���̍��W

				const double temp_x_coil_1 = output_coil[0].coil_pos_x * output_coil[0].COS_XY + output_coil[0].coil_pos_y * output_coil[0].SIN_XY;		//���W�ϊ���̍��W
				const double temp_y_coil_1 = -output_coil[0].coil_pos_x * output_coil[0].SIN_XY + output_coil[0].coil_pos_y * output_coil[0].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_coil_1 = output_coil[0].coil_pos_z;																					//���W�ϊ���̍��W
				const double temp_x_coil_2 = temp_x_coil_1 * output_coil[0].COS_XZ + temp_z_coil_1 * output_coil[0].SIN_XZ;		//���W�ϊ���̍��W
				const double temp_y_coil_2 = temp_y_coil_1;																		//���W�ϊ���̍��W
				const double temp_z_coil_2 = -temp_x_coil_1 * output_coil[0].SIN_XZ + temp_z_coil_1 * output_coil[0].COS_XZ;	//���W�ϊ���̍��W


				const double temp_z_rel = temp_x_3 - (temp_x_coil_2 - endpoint);						//���W�ϊ���̃R�C����[����̋���
				const double temp_r_rel = ((temp_y_3 - temp_y_coil_2) * (temp_y_3 - temp_y_coil_2)) 
										+ ((temp_z_3 - temp_z_coil_2) * (temp_z_3 - temp_z_coil_2));    //���W�ϊ���̃R�C�����S������̋�����2��

				const int z_num = (temp_z_rel == 0) ? 0 : (int)(temp_z_rel / ElmSize + 0.5);		//�~�����W�n�ɂ�������W
				const int r_num = (temp_r_rel == 0) ? 0 : (int)(sqrt(temp_r_rel) / ElmSize + 0.5);	//�~�����W�n�ɂ�������W
				const int n = z_num * WidthNum + r_num;

				if ((z_num < 0) || (z_num >= HeightNum) || (r_num >= WidthNum)){
					Output_B_x_0[m] = 0.0;
					Output_B_y_0[m] = 0.0;
					Output_B_z_0[m] = 0.0;
				}
				else{
					const double radian = atan2((temp_z_3 - temp_z_coil_2), (temp_y_3 - temp_y_coil_2));
					const double b_x = B_z[n];
					const double b_y = cos(radian) * B_r[n];
					const double b_z = sin(radian) * B_r[n];
					const double output_x = b_x * output_coil[0].COS_XZ - b_z * output_coil[0].SIN_XZ;
					const double output_y = b_y;
					const double output_z = b_x * output_coil[0].SIN_XZ + b_z * output_coil[0].COS_XZ;
					Output_B_x_0[m] = output_x * output_coil[0].COS_XY - output_y * output_coil[0].SIN_XY;
					Output_B_y_0[m] = output_x * output_coil[0].SIN_XY + output_y * output_coil[0].COS_XY;
					Output_B_z_0[m] = output_z;
				}
			}
		}
	}


	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				const int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				const double temp_x_1 = i * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_y_1 = j * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_z_1 = k * ElmSize_output;		//���W�ϊ���̍��W

				const double temp_x_2 = temp_x_1 * output_coil[1].COS_XY + temp_y_1 * output_coil[1].SIN_XY;	//���W�ϊ���̍��W
				const double temp_y_2 = -temp_x_1 * output_coil[1].SIN_XY + temp_y_1 * output_coil[1].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_2 = temp_z_1;																//���W�ϊ���̍��W
				const double temp_x_3 = temp_x_2 * output_coil[1].COS_XZ + temp_z_2 * output_coil[1].SIN_XZ;	//���W�ϊ���̍��W
				const double temp_y_3 = temp_y_2;																//���W�ϊ���̍��W
				const double temp_z_3 = -temp_x_2 * output_coil[1].SIN_XZ + temp_z_2 * output_coil[1].COS_XZ;	//���W�ϊ���̍��W

				const double temp_x_coil_1 = output_coil[1].coil_pos_x * output_coil[1].COS_XY + output_coil[1].coil_pos_y * output_coil[1].SIN_XY;		//���W�ϊ���̍��W
				const double temp_y_coil_1 = -output_coil[1].coil_pos_x * output_coil[1].SIN_XY + output_coil[1].coil_pos_y * output_coil[1].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_coil_1 = output_coil[1].coil_pos_z;																					//���W�ϊ���̍��W
				const double temp_x_coil_2 = temp_x_coil_1 * output_coil[1].COS_XZ + temp_z_coil_1 * output_coil[1].SIN_XZ;		//���W�ϊ���̍��W
				const double temp_y_coil_2 = temp_y_coil_1;																		//���W�ϊ���̍��W
				const double temp_z_coil_2 = -temp_x_coil_1 * output_coil[1].SIN_XZ + temp_z_coil_1 * output_coil[1].COS_XZ;	//���W�ϊ���̍��W


				const double temp_z_rel = temp_x_3 - (temp_x_coil_2 - endpoint);						//���W�ϊ���̃R�C����[����̋���
				const double temp_r_rel = ((temp_y_3 - temp_y_coil_2) * (temp_y_3 - temp_y_coil_2))
					+ ((temp_z_3 - temp_z_coil_2) * (temp_z_3 - temp_z_coil_2));    //���W�ϊ���̃R�C�����S������̋�����2��

				const int z_num = (temp_z_rel == 0) ? 0 : (int)(temp_z_rel / ElmSize + 0.5);		//�~�����W�n�ɂ�������W
				const int r_num = (temp_r_rel == 0) ? 0 : (int)(sqrt(temp_r_rel) / ElmSize + 0.5);	//�~�����W�n�ɂ�������W
				const int n = z_num * WidthNum + r_num;

				if ((z_num < 0) || (z_num >= HeightNum) || (r_num >= WidthNum)){
					Output_B_x_1[m] = 0.0;
					Output_B_y_1[m] = 0.0;
					Output_B_z_1[m] = 0.0;
				}
				else{
					const double radian = atan2((temp_z_3 - temp_z_coil_2), (temp_y_3 - temp_y_coil_2));
					const double b_x = B_z[n];
					const double b_y = cos(radian) * B_r[n];
					const double b_z = sin(radian) * B_r[n];
					const double output_x = b_x * output_coil[1].COS_XZ - b_z * output_coil[1].SIN_XZ;
					const double output_y = b_y;
					const double output_z = b_x * output_coil[1].SIN_XZ + b_z * output_coil[1].COS_XZ;
					Output_B_x_1[m] = output_x * output_coil[1].COS_XY - output_y * output_coil[1].SIN_XY;
					Output_B_y_1[m] = output_x * output_coil[1].SIN_XY + output_y * output_coil[1].COS_XY;
					Output_B_z_1[m] = output_z;
				}
			}
		}
	}

	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				const int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				const double temp_x_1 = i * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_y_1 = j * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_z_1 = k * ElmSize_output;		//���W�ϊ���̍��W

				const double temp_x_2 = temp_x_1 * output_coil[2].COS_XY + temp_y_1 * output_coil[2].SIN_XY;	//���W�ϊ���̍��W
				const double temp_y_2 = -temp_x_1 * output_coil[2].SIN_XY + temp_y_1 * output_coil[2].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_2 = temp_z_1;																//���W�ϊ���̍��W
				const double temp_x_3 = temp_x_2 * output_coil[2].COS_XZ + temp_z_2 * output_coil[2].SIN_XZ;	//���W�ϊ���̍��W
				const double temp_y_3 = temp_y_2;																//���W�ϊ���̍��W
				const double temp_z_3 = -temp_x_2 * output_coil[2].SIN_XZ + temp_z_2 * output_coil[2].COS_XZ;	//���W�ϊ���̍��W

				const double temp_x_coil_1 = output_coil[2].coil_pos_x * output_coil[2].COS_XY + output_coil[2].coil_pos_y * output_coil[2].SIN_XY;		//���W�ϊ���̍��W
				const double temp_y_coil_1 = -output_coil[2].coil_pos_x * output_coil[2].SIN_XY + output_coil[2].coil_pos_y * output_coil[2].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_coil_1 = output_coil[2].coil_pos_z;																					//���W�ϊ���̍��W
				const double temp_x_coil_2 = temp_x_coil_1 * output_coil[2].COS_XZ + temp_z_coil_1 * output_coil[2].SIN_XZ;		//���W�ϊ���̍��W
				const double temp_y_coil_2 = temp_y_coil_1;																		//���W�ϊ���̍��W
				const double temp_z_coil_2 = -temp_x_coil_1 * output_coil[2].SIN_XZ + temp_z_coil_1 * output_coil[2].COS_XZ;	//���W�ϊ���̍��W


				const double temp_z_rel = temp_x_3 - (temp_x_coil_2 - endpoint);						//���W�ϊ���̃R�C����[����̋���
				const double temp_r_rel = ((temp_y_3 - temp_y_coil_2) * (temp_y_3 - temp_y_coil_2))
					+ ((temp_z_3 - temp_z_coil_2) * (temp_z_3 - temp_z_coil_2));    //���W�ϊ���̃R�C�����S������̋�����2��

				const int z_num = (temp_z_rel == 0) ? 0 : (int)(temp_z_rel / ElmSize + 0.5);		//�~�����W�n�ɂ�������W
				const int r_num = (temp_r_rel == 0) ? 0 : (int)(sqrt(temp_r_rel) / ElmSize + 0.5);	//�~�����W�n�ɂ�������W
				const int n = z_num * WidthNum + r_num;

				if ((z_num < 0) || (z_num >= HeightNum) || (r_num >= WidthNum)){
					Output_B_x_2[m] = 0.0;
					Output_B_y_2[m] = 0.0;
					Output_B_z_2[m] = 0.0;
				}
				else{
					const double radian = atan2((temp_z_3 - temp_z_coil_2), (temp_y_3 - temp_y_coil_2));
					const double b_x = B_z[n];
					const double b_y = cos(radian) * B_r[n];
					const double b_z = sin(radian) * B_r[n];
					const double output_x = b_x * output_coil[2].COS_XZ - b_z * output_coil[2].SIN_XZ;
					const double output_y = b_y;
					const double output_z = b_x * output_coil[2].SIN_XZ + b_z * output_coil[2].COS_XZ;
					Output_B_x_2[m] = output_x * output_coil[2].COS_XY - output_y * output_coil[2].SIN_XY;
					Output_B_y_2[m] = output_x * output_coil[2].SIN_XY + output_y * output_coil[2].COS_XY;
					Output_B_z_2[m] = output_z;
				}
			}
		}
	}

	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				const int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				const double temp_x_1 = i * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_y_1 = j * ElmSize_output;		//���W�ϊ��O�̍��W
				const double temp_z_1 = k * ElmSize_output;		//���W�ϊ���̍��W

				const double temp_x_2 = temp_x_1 * output_coil[3].COS_XY + temp_y_1 * output_coil[3].SIN_XY;	//���W�ϊ���̍��W
				const double temp_y_2 = -temp_x_1 * output_coil[3].SIN_XY + temp_y_1 * output_coil[3].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_2 = temp_z_1;																//���W�ϊ���̍��W
				const double temp_x_3 = temp_x_2 * output_coil[3].COS_XZ + temp_z_2 * output_coil[3].SIN_XZ;	//���W�ϊ���̍��W
				const double temp_y_3 = temp_y_2;																//���W�ϊ���̍��W
				const double temp_z_3 = -temp_x_2 * output_coil[3].SIN_XZ + temp_z_2 * output_coil[3].COS_XZ;	//���W�ϊ���̍��W

				const double temp_x_coil_1 = output_coil[3].coil_pos_x * output_coil[3].COS_XY + output_coil[3].coil_pos_y * output_coil[3].SIN_XY;		//���W�ϊ���̍��W
				const double temp_y_coil_1 = -output_coil[3].coil_pos_x * output_coil[3].SIN_XY + output_coil[3].coil_pos_y * output_coil[3].COS_XY;	//���W�ϊ���̍��W
				const double temp_z_coil_1 = output_coil[3].coil_pos_z;																					//���W�ϊ���̍��W
				const double temp_x_coil_2 = temp_x_coil_1 * output_coil[3].COS_XZ + temp_z_coil_1 * output_coil[3].SIN_XZ;		//���W�ϊ���̍��W
				const double temp_y_coil_2 = temp_y_coil_1;																		//���W�ϊ���̍��W
				const double temp_z_coil_2 = -temp_x_coil_1 * output_coil[3].SIN_XZ + temp_z_coil_1 * output_coil[3].COS_XZ;	//���W�ϊ���̍��W


				const double temp_z_rel = temp_x_3 - (temp_x_coil_2 - endpoint);						//���W�ϊ���̃R�C����[����̋���
				const double temp_r_rel = ((temp_y_3 - temp_y_coil_2) * (temp_y_3 - temp_y_coil_2))
					+ ((temp_z_3 - temp_z_coil_2) * (temp_z_3 - temp_z_coil_2));    //���W�ϊ���̃R�C�����S������̋�����2��

				const int z_num = (temp_z_rel == 0) ? 0 : (int)(temp_z_rel / ElmSize + 0.5);		//�~�����W�n�ɂ�������W
				const int r_num = (temp_r_rel == 0) ? 0 : (int)(sqrt(temp_r_rel) / ElmSize + 0.5);	//�~�����W�n�ɂ�������W
				const int n = z_num * WidthNum + r_num;

				if ((z_num < 0) || (z_num >= HeightNum) || (r_num >= WidthNum)){
					Output_B_x_3[m] = 0.0;
					Output_B_y_3[m] = 0.0;
					Output_B_z_3[m] = 0.0;
				}
				else{
					const double radian = atan2((temp_z_3 - temp_z_coil_2), (temp_y_3 - temp_y_coil_2));
					const double b_x = B_z[n];
					const double b_y = cos(radian) * B_r[n];
					const double b_z = sin(radian) * B_r[n];
					const double output_x = b_x * output_coil[3].COS_XZ - b_z * output_coil[3].SIN_XZ;
					const double output_y = b_y;
					const double output_z = b_x * output_coil[3].SIN_XZ + b_z * output_coil[3].COS_XZ;
					Output_B_x_3[m] = output_x * output_coil[3].COS_XY - output_y * output_coil[3].SIN_XY;
					Output_B_y_3[m] = output_x * output_coil[3].SIN_XY + output_y * output_coil[3].COS_XY;
					Output_B_z_3[m] = output_z;
				}
			}
		}
	}

}

void Magnetic_Calc::Cal_Temp(char *filename_0, char *filename_1, char *filename_2, char *filename_3){

	ofstream fout_B0;
	ofstream fout_B1;
	ofstream fout_B2;
	ofstream fout_B3;

	fout_B0.open(filename_0, ios_base::out | ios_base::trunc);
	fout_B1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_B2.open(filename_2, ios_base::out | ios_base::trunc);
	fout_B3.open(filename_3, ios_base::out | ios_base::trunc);

	fout_B0	<< "Pos_x [m],Pos_y [m],Pos_z [m],B_x [T],B_y [T],B_z [T]"<< endl;

	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				fout_B0 << i * ElmSize_output << "," << j * ElmSize_output << "," << k * ElmSize_output << "," << Output_B_x_0[m] << "," << Output_B_y_0[m] << "," << Output_B_z_0[m] << '\n';
				fout_B1 << i * ElmSize_output << "," << j * ElmSize_output << "," << k * ElmSize_output << "," << Output_B_x_1[m] << "," << Output_B_y_1[m] << "," << Output_B_z_1[m] << '\n';
				fout_B2 << i * ElmSize_output << "," << j * ElmSize_output << "," << k * ElmSize_output << "," << Output_B_x_2[m] << "," << Output_B_y_2[m] << "," << Output_B_z_2[m] << '\n';
				fout_B3 << i * ElmSize_output << "," << j * ElmSize_output << "," << k * ElmSize_output << "," << Output_B_x_3[m] << "," << Output_B_y_3[m] << "," << Output_B_z_3[m] << '\n';
			}
		}
	}

	fout_B0.close();
	fout_B1.close();
	fout_B2.close();
	fout_B3.close();

}


void Magnetic_Calc::OutPut_Coil(char *filename_1, char *filename_2){

	cout << "�R�C���ʒu���o�͂��܂�" << endl;

	ofstream fout_potential1;
	ofstream fout_potential2;

	fout_potential1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_potential2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_potential1
		<< "# AVS field file\n" << "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum + i;
			fout_potential2 << CurrentDensity[m] << '\n';
		}
	}

	fout_potential1.close();
	fout_potential2.close();
}

void Magnetic_Calc::OutPut_Potential(char *filename_1, char *filename_2){

	cout << "�|�e���V�����f�[�^���o�͂��܂�" << endl;

	ofstream fout_potential1;
	ofstream fout_potential2;

	fout_potential1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_potential2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_potential1
		<< "# AVS field file\n" << "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";
	
	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum + i;
			fout_potential2 << Potential[m] << '\n';
		}
	}

	fout_potential1.close();
	fout_potential2.close();
}

void Magnetic_Calc::OutPut_B_1(char *filename_1, char *filename_2){

	cout << "�������x�f�[�^���o�͂��܂�" << endl;

	ofstream fout_B1;
	ofstream fout_B2;

	fout_B1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_B2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_B1
		<< "# AVS field file\n" << "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";
	
	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){
			int m = k * WidthNum + i;
			fout_B2 << sqrt(pow(B_r[m], 2) + pow(B_z[m], 2)) << '\n';
		}
	}

	fout_B1.close();
	fout_B2.close();
}

void Magnetic_Calc::OutPut_B_2(char *filename_1, char *filename_2){

	cout << "�������x�f�[�^���o�͂��܂� ������I" << endl;

	ofstream fout_B1;
	ofstream fout_B2;

	fout_B1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_B2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_B1
		<< "# AVS field file\n" << "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n" << "veclen = 2\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii" << ' ' << "offset=0" << ' ' << "stride=2\n"
		<< "variable 2 file=" << filename_2 << ' ' << "filetype=ascii" << ' ' << "offset=1" << ' ' << "stride=2\n";
	
	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){
			int m = k * WidthNum + i;
			fout_B2 << B_r[m] << " " <<B_z[m] << '\n';
		}
	}

	fout_B1.close();
	fout_B2.close();
}

void Magnetic_Calc::OutPut_B(char *filename_1, char *filename_2){

	cout << "�������x�f�[�^���o�͂��܂�" << endl;

	ofstream fout_B_1;
	ofstream fout_B_2;

	fout_B_1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_B_2.open(filename_2, ios_base::out | ios_base::trunc);

	for (int i = 0; i<WidthNum; i++){
		int m = (HeightNum / 2) * WidthNum + i;
		fout_B_1 << i * ElmSize << "," << B_z[m] << '\n';
	}

	for (int k = 0; k<HeightNum; k++){
		int m = k * WidthNum + 0;
		fout_B_2 << k * ElmSize << "," << B_z[m] << '\n';
	}

	fout_B_1.close();
	fout_B_2.close();
}

void Magnetic_Calc::OutPut_B_Cylindrical(char *filename_0_a, char *filename_0_b,
	char *filename_1_a, char *filename_1_b,
	char *filename_2_a, char *filename_2_b,
	char *filename_3_a, char *filename_3_b){

	cout << "�f�J���g���W�n�̎������x���o��" << endl;

	ofstream fout_B0_a;
	ofstream fout_B0_b;
	ofstream fout_B1_a;
	ofstream fout_B1_b;
	ofstream fout_B2_a;
	ofstream fout_B2_b;
	ofstream fout_B3_a;
	ofstream fout_B3_b;

	fout_B0_a.open(filename_0_a, ios_base::out | ios_base::trunc);
	fout_B0_b.open(filename_0_b, ios_base::out | ios_base::trunc);
	fout_B1_a.open(filename_1_a, ios_base::out | ios_base::trunc);
	fout_B1_b.open(filename_1_b, ios_base::out | ios_base::trunc);
	fout_B2_a.open(filename_2_a, ios_base::out | ios_base::trunc);
	fout_B2_b.open(filename_2_b, ios_base::out | ios_base::trunc);
	fout_B3_a.open(filename_3_a, ios_base::out | ios_base::trunc);
	fout_B3_b.open(filename_3_b, ios_base::out | ios_base::trunc);

	fout_B0_a
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum_output << '\n'
		<< "dim2 = " << DepthNum_output << '\n'
		<< "dim3 = " << HeightNum_output << '\n'
		<< "nspace = 3\n" << "veclen = 3\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_0_b << ' ' << "filetype=ascii" << ' ' << "offset=0" << ' ' << "stride=3\n"
		<< "variable 2 file=" << filename_0_b << ' ' << "filetype=ascii" << ' ' << "offset=1" << ' ' << "stride=3\n"
		<< "variable 3 file=" << filename_0_b << ' ' << "filetype=ascii" << ' ' << "offset=2" << ' ' << "stride=3\n";

	fout_B1_a
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum_output << '\n'
		<< "dim2 = " << DepthNum_output << '\n'
		<< "dim3 = " << HeightNum_output << '\n'
		<< "nspace = 3\n" << "veclen = 3\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_1_b << ' ' << "filetype=ascii" << ' ' << "offset=0" << ' ' << "stride=3\n"
		<< "variable 2 file=" << filename_1_b << ' ' << "filetype=ascii" << ' ' << "offset=1" << ' ' << "stride=3\n"
		<< "variable 3 file=" << filename_1_b << ' ' << "filetype=ascii" << ' ' << "offset=2" << ' ' << "stride=3\n";

	fout_B2_a
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum_output << '\n'
		<< "dim2 = " << DepthNum_output << '\n'
		<< "dim3 = " << HeightNum_output << '\n'
		<< "nspace = 3\n" << "veclen = 3\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2_b << ' ' << "filetype=ascii" << ' ' << "offset=0" << ' ' << "stride=3\n"
		<< "variable 2 file=" << filename_2_b << ' ' << "filetype=ascii" << ' ' << "offset=1" << ' ' << "stride=3\n"
		<< "variable 3 file=" << filename_2_b << ' ' << "filetype=ascii" << ' ' << "offset=2" << ' ' << "stride=3\n";

	fout_B3_a
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum_output << '\n'
		<< "dim2 = " << DepthNum_output << '\n'
		<< "dim3 = " << HeightNum_output << '\n'
		<< "nspace = 3\n" << "veclen = 3\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_3_b << ' ' << "filetype=ascii" << ' ' << "offset=0" << ' ' << "stride=3\n"
		<< "variable 2 file=" << filename_3_b << ' ' << "filetype=ascii" << ' ' << "offset=1" << ' ' << "stride=3\n"
		<< "variable 3 file=" << filename_3_b << ' ' << "filetype=ascii" << ' ' << "offset=2" << ' ' << "stride=3\n";


	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				fout_B0_b << Output_B_x_0[m] << " " << Output_B_y_0[m] << " " << Output_B_z_0[m] << '\n';
				fout_B1_b << Output_B_x_1[m] << " " << Output_B_y_1[m] << " " << Output_B_z_1[m] << '\n';
				fout_B2_b << Output_B_x_2[m] << " " << Output_B_y_2[m] << " " << Output_B_z_2[m] << '\n';
				fout_B3_b << Output_B_x_3[m] << " " << Output_B_y_3[m] << " " << Output_B_z_3[m] << '\n';
			}
		}
	}

	fout_B0_a.close();
	fout_B0_b.close();	
	fout_B1_a.close();
	fout_B1_b.close();
	fout_B2_a.close();
	fout_B2_b.close();
	fout_B3_a.close();
	fout_B3_b.close();
}

void Magnetic_Calc::OutPut_MagFieldBinary(char *filename_0, char *filename_1, char *filename_2, char *filename_3){

	cout << "�����v�Z�p�̎��E�f�[�^���o�͂��܂�" << endl;

	ofstream ofs_0;
	ofs_0.open(filename_0, ios_base::out | ios_base::binary | ios_base::trunc);
	ofstream ofs_1;
	ofs_1.open(filename_1, ios_base::out | ios_base::binary | ios_base::trunc);
	ofstream ofs_2;
	ofs_2.open(filename_2, ios_base::out | ios_base::binary | ios_base::trunc);
	ofstream ofs_3;
	ofs_3.open(filename_3, ios_base::out | ios_base::binary | ios_base::trunc);

	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				ofs_0.write((const char*)&Output_B_x_0[m], sizeof(double));
				ofs_1.write((const char*)&Output_B_x_1[m], sizeof(double));
				ofs_2.write((const char*)&Output_B_x_2[m], sizeof(double));
				ofs_3.write((const char*)&Output_B_x_3[m], sizeof(double));
			}
		}
	}

	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				ofs_0.write((const char*)&Output_B_y_0[m], sizeof(double));
				ofs_1.write((const char*)&Output_B_y_1[m], sizeof(double));
				ofs_2.write((const char*)&Output_B_y_2[m], sizeof(double));
				ofs_3.write((const char*)&Output_B_y_3[m], sizeof(double));

			}
		}
	}

	for (int k = 0; k<HeightNum_output; k++){
		for (int j = 0; j<DepthNum_output; j++){
			for (int i = 0; i < WidthNum_output; i++){

				int m = k * WidthNum_output * DepthNum_output + j * WidthNum_output + i;
				ofs_0.write((const char*)&Output_B_z_0[m], sizeof(double));
				ofs_1.write((const char*)&Output_B_z_1[m], sizeof(double));
				ofs_2.write((const char*)&Output_B_z_2[m], sizeof(double));
				ofs_3.write((const char*)&Output_B_z_3[m], sizeof(double));
			}
		}
	}

	ofs_0.close();
	ofs_1.close();
	ofs_2.close();
	ofs_3.close();

}

void Magnetic_Calc::OutPut_Condition(char *filename){

	ofstream ofs;
	ofs.open(filename, ios_base::out | ios_base::trunc);
	
	ofs << "mesh_size =  " << ElmSize << " [m] " << endl;
	ofs << "areaWidth =  " << areaWidth << " [m] " << endl;
	ofs << "areaHeight =  " << areaHeight << " [m] " << endl;
	ofs << "WidthNum =  " << WidthNum << endl;
	ofs << "HeightNum =  " << HeightNum << endl;
	ofs << "nodeNum =  " << nodeNum << endl;
	ofs << "endpoint =  " << endpoint << " [m] " << endl;
	ofs << "wire_diameter =  " << wire_diameter << " [m] " << endl;
	ofs << "wire_length =  " << wire_totallength << " [m] " << endl;
	ofs << "coil_length =  " << coil_length << " [m] " << endl;
	ofs << "coil_diameter =  " << coil_diameter << " [m] " << endl;
	ofs << "wire_turn =  " << wire_turn << endl;
	ofs << "turn_length =  " << turn_number_height << endl;
	ofs << "permeability =  " << permeability << " [H/m] " << endl;
	ofs << "I =  " << I << " [A] " << endl;
	ofs << "i =  " << I / ((wire_diameter * 0.5) * (wire_diameter * 0.5) * PI ) << " [A/m2] " << endl;
	ofs << "p =  " << p << " [��m] " << endl;
	ofs << "R =  " << R << " [��] " << endl;
	ofs << "W =  " << W << " [J/s] " << endl;

	ofs << endl << "for particle simulation as below" << endl;

	ofs << "mesh_size_sim =  " << ElmSize_output << " [m] " << endl;
	ofs << "areaWidth_sim =  " << areaWidth_output << " [m] " << endl;
	ofs << "areaDepth_sim =  " << areaDepth_output << " [m] " << endl;
	ofs << "areaHeight_sim =  " << areaHeight_output << " [m] " << endl;
	ofs << "WidthNum_sim =  " << WidthNum_output << endl;
	ofs << "DepthNum_sim =  " << DepthNum_output << endl;
	ofs << "HeightNum_sim =  " << HeightNum_output << endl;
	ofs << "nodeNum_sim =  " << nodeNum_output << endl;

	ofs << "incline_coil_XY_0 =  " << output_coil[0].incline_XY / (PI / 180) << " [��] " << endl;
	ofs << "incline_coil_XZ_0 =  " << output_coil[0].incline_XZ / (PI / 180) << " [��] " << endl;
	ofs << "coil_pos_X_0 =  " << output_coil[0].coil_pos_x << " [m] " << endl;
	ofs << "coil_pos_Y_0 =  " << output_coil[0].coil_pos_y << " [m] " << endl;
	ofs << "coil_pos_Z_0 =  " << output_coil[0].coil_pos_z << " [m] " << endl;

	ofs << "incline_coil_XY_1 =  " << output_coil[1].incline_XY / (PI / 180) << " [��] " << endl;
	ofs << "incline_coil_XZ_1 =  " << output_coil[1].incline_XZ / (PI / 180) << " [��] " << endl;
	ofs << "coil_pos_X_1 =  " << output_coil[1].coil_pos_x << " [m] " << endl;
	ofs << "coil_pos_Y_1 =  " << output_coil[1].coil_pos_y << " [m] " << endl;
	ofs << "coil_pos_Z_1 =  " << output_coil[1].coil_pos_z << " [m] " << endl;

	ofs << "incline_coil_XY_2 =  " << output_coil[2].incline_XY / (PI / 180) << " [��] " << endl;
	ofs << "incline_coil_XZ_2 =  " << output_coil[2].incline_XZ / (PI / 180) << " [��] " << endl;
	ofs << "coil_pos_X_2 =  " << output_coil[2].coil_pos_x << " [m] " << endl;
	ofs << "coil_pos_Y_2 =  " << output_coil[2].coil_pos_y << " [m] " << endl;
	ofs << "coil_pos_Z_2 =  " << output_coil[2].coil_pos_z << " [m] " << endl;

	ofs << "incline_coil_XY_3 =  " << output_coil[3].incline_XY / (PI / 180) << " [��] " << endl;
	ofs << "incline_coil_XZ_3 =  " << output_coil[3].incline_XZ / (PI / 180) << " [��] " << endl;
	ofs << "coil_pos_X_3 =  " << output_coil[3].coil_pos_x << " [m] " << endl;
	ofs << "coil_pos_Y_3 =  " << output_coil[3].coil_pos_y << " [m] " << endl;
	ofs << "coil_pos_Z_3 =  " << output_coil[3].coil_pos_z << " [m] " << endl;

	ofs.close();
}