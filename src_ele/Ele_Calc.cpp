#include "Ele_Calc.h"

Electric_Calc::Electric_Calc() :PI(3.14159265358979323846264338328){		//�R���X�g���N�^�D�ϐ��̏����������Ă���D

	//��{�p�����[�^
	ElmSize = 4.0e-6;			//���b�V���T�C�Y
	areaWidth = 5200.0e-6;		//�w�����v�Z�̈�
	areaHeight = 12032.0e-6;	//�y�����v�Z�̈�	
	Type_boundary_X = 1;		//���̋��E�����̐ݒ�@0  = �m�C�}�� / 1 = ���� ����I��
	
	spElePerm = 8.854187816e-12;		//�U�d��
	AirRelElePerm = 1.00059;			//��Ԃ̔�U�d��
	Volt = 750;							//����d��		
	TFC = 1.0e-3;						//��������
	SOR = 1.8;							//�����W��	
	OMP_CALC = true;					//���񉻂��邩
	thread = 6;							//����

	//���d�ɂ�����ꍇ
	cylinder_Num = 0;					//�d�ɐ�
	CylinderElePerm = 5.0;				//�U�d��

	if (cylinder_Num){
		cylinder.resize(cylinder_Num);		

		for (int i = 0; i < cylinder_Num; i++){

			cylinder[i].Pos_X = 500.0e-6 + 1000.0e-6 *i;
			cylinder[i].Pos_Y = 0;
			cylinder[i].Pos_Z = 7000.0e-6;
			cylinder[i].pos_X = (int)(cylinder[i].Pos_X / ElmSize + 0.5);
			cylinder[i].pos_Y = (int)(cylinder[i].Pos_Y / ElmSize + 0.5);
			cylinder[i].pos_Z = (int)(cylinder[i].Pos_Z / ElmSize + 0.5);

			cylinder[i].Radius = 100.0e-6;				  //�d�ɔ��a
			cylinder[i].Dielectric_Radius = 200.0e-6;     //�d�ɔ핞���܂߂����a
			cylinder[i].radius = (int)(cylinder[i].Radius / ElmSize + 0.5);
			cylinder[i].diele_radius = (int)(cylinder[i].Dielectric_Radius / ElmSize + 0.5);
		}
	}

	//�����̂�����ꍇ
	wallNum_electrode = 4;

	if (wallNum_electrode){
		wall_electrode.resize(wallNum_electrode);

		for (int i = 0; i < wallNum_electrode; i++){

			wall_electrode[i].Pos_X[0] = 500.0e-6 + 1300.0e-6 * i;
			wall_electrode[i].Pos_Y[0] = 0;
			wall_electrode[i].Pos_Z[0] = 6000.0e-6;

			wall_electrode[i].Pos_X[1] = 800.0e-6 + 1300.0e-6 * i;
			wall_electrode[i].Pos_Y[1] = 0;
			wall_electrode[i].Pos_Z[1] = 6000.0e-6;

			wall_electrode[i].Pos_X[2] = 300.0e-6 + 1300.0e-6 * i;
			wall_electrode[i].Pos_Y[2] = 0;
			wall_electrode[i].Pos_Z[2] = 6020.0e-6;

			wall_electrode[i].Pos_X[3] = 800.0e-6 + 1300.0e-6 * i;
			wall_electrode[i].Pos_Y[3] = 0;
			wall_electrode[i].Pos_Z[3] = 6020.0e-6;

			for (int n = 0; n < 4; n++){
				wall_electrode[i].pos_X[n] = (int)(wall_electrode[i].Pos_X[n] / ElmSize + 0.5);
				wall_electrode[i].pos_Y[n] = (int)(wall_electrode[i].Pos_Y[n] / ElmSize + 0.5);
				wall_electrode[i].pos_Z[n] = (int)(wall_electrode[i].Pos_Z[n] / ElmSize + 0.5);
			}
		}
	}

	//�����̂�����ꍇ
	wallNum_dielectric = 1;	
	wallElePerm = 3.1;

	if (wallNum_dielectric){
		wall_dielectric.resize(wallNum_dielectric);

		for (int i = 0; i < wallNum_dielectric; i++){

			wall_dielectric[i].Pos_X[0] = 0;
			wall_dielectric[i].Pos_Y[0] = 0;
			wall_dielectric[i].Pos_Z[0] = 0;

			wall_dielectric[i].Pos_X[1] = areaWidth;
			wall_dielectric[i].Pos_Y[1] = 0;
			wall_dielectric[i].Pos_Z[1] = 0;

			wall_dielectric[i].Pos_X[2] = 0;
			wall_dielectric[i].Pos_Y[2] = 0;
			wall_dielectric[i].Pos_Z[2] = 6132.0e-6;

			wall_dielectric[i].Pos_X[3] = areaWidth;
			wall_dielectric[i].Pos_Y[3] = 0;
			wall_dielectric[i].Pos_Z[3] = 6132.0e-6;

			for (int n = 0; n < 4; n++){
				wall_dielectric[i].pos_X[n] = (int)(wall_dielectric[i].Pos_X[n] / ElmSize + 0.5);
				wall_dielectric[i].pos_Y[n] = (int)(wall_dielectric[i].Pos_Y[n] / ElmSize + 0.5);
				wall_dielectric[i].pos_Z[n] = (int)(wall_dielectric[i].Pos_Z[n] / ElmSize + 0.5);
			}
		}
	}
	
	//�ȍ~�͕ϐ��̏������E�p�����[�^�̗��U��
	WidthNum = (int)(areaWidth / ElmSize + 0.5);
	HeightNum = (int)(areaHeight / ElmSize + 0.5);
	cellNum = WidthNum * HeightNum;
	if (cellNum == 0)	cellNum = 1;
	convergentCount = 0;

	ElePerm.assign(cellNum, 0.0);
	Potential.assign(cellNum, 0.0);
	E_x.assign(cellNum, 0.0);
	E_z.assign(cellNum, 0.0);
	bcc.assign(cellNum, 0);
	fcc.assign(cellNum, 0);
	ccb.assign(cellNum, 0);
	ccf.assign(cellNum, 0);
	bcb.assign(cellNum, 0);
	fcb.assign(cellNum, 0);
	bcf.assign(cellNum, 0);
	fcf.assign(cellNum, 0);
	CalcFlag = new bool[cellNum];
	for (int i = 0; i < cellNum; i++)	CalcFlag[i] = true;

	cout << "�v�Z������ݒ芮��" << endl;
}

Electric_Calc::~Electric_Calc(){

	delete[] CalcFlag;
}

void Electric_Calc::Setting_Boundary_Condition_V(){

	//���̓d�ʂɊւ��鋫�E�����ݒ�
	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			bcc[m] = m - 1;
			fcc[m] = m + 1;
			ccb[m] = m - WidthNum;
			ccf[m] = m + WidthNum;
		}
	}

	if (Type_boundary_X == 0){			Setting_Boundary_Condition_V_X_n();
	}
	else if (Type_boundary_X == 1){		Setting_Boundary_Condition_V_X_r();
	}

	Setting_Boundary_Condition_V_Z_n();
	
	cout << "�d�ʂɊւ���[�̋��E�����ݒ�" << endl;
}

void Electric_Calc::Setting_Boundary_Condition_V_X_n(){
	
	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			if (i == 0)					bcc[m] = m + 1;
			if (i == WidthNum - 1)		fcc[m] = m - 1;
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_X_r(){

	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			if (i == 0)					bcc[m] = m + WidthNum - 1;
			if (i == WidthNum - 1)		fcc[m] = m - WidthNum + 1;
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			if (k == 0)					ccb[m] = m + WidthNum;
			if (k == HeightNum - 1)		ccf[m] = m - WidthNum;
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e(){

	//���̗U�d���Ɋւ��鋫�E�����ݒ�
	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;
			bcb[m] = m - WidthNum - 1;
			fcb[m] = m - WidthNum;
			bcf[m] = m - 1;
			fcf[m] = m;
		}
	}

	if (Type_boundary_X == 0)			Setting_Boundary_Condition_e_X_n_Z_n();
	else if (Type_boundary_X == 1)		Setting_Boundary_Condition_e_X_r_Z_n();
	
	cout << "�U�d���Ɋւ���[�̋��E�����ݒ�" << endl;
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){
			int m = k * WidthNum + i;

			////////////////////////////////////////////////////////
			if (i == 0){
				bcb[m] = m - WidthNum;
				bcf[m] = m;
			}
			if (i == WidthNum - 1){
				fcb[m] = m - WidthNum - 1;
				fcf[m] = m - 1;
			}
			if (k == 0){
				bcb[m] = m - 1;
				fcb[m] = m;
			}
			if (k == HeightNum - 1){
				bcf[m] = m - WidthNum - 1;
				fcf[m] = m - WidthNum;
			}
			/////////////////////////////////////////////////////////
			if ((i == 0) && (k == 0)){
				bcb[m] = m;
				fcb[m] = m;
				bcf[m] = m;
			}
			if ((i == WidthNum - 1) && (k == 0)){
				bcb[m] = m - 1;
				fcb[m] = m - 1;
				fcf[m] = m - 1;
			}
			if ((i == 0) && (k == HeightNum - 1)){
				bcb[m] = m - WidthNum;
				bcf[m] = m - WidthNum;
				fcf[m] = m - WidthNum;
			}
			if ((i == WidthNum - 1) && (k == HeightNum - 1)){
				fcb[m] = m - WidthNum - 1;
				bcf[m] = m - WidthNum - 1;
				fcf[m] = m - WidthNum - 1;
			}
			/////////////////////////////////////////////////////////
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_X_r_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){
			int m = k * WidthNum + i;

			////////////////////////////////////////////////////////
			if (i == 0){
				bcb[m] = m - 2;
				bcf[m] = m + WidthNum - 2;
			}
			if (i == WidthNum - 1){
				fcb[m] = m - WidthNum - WidthNum + 1;
				fcf[m] = m - WidthNum + 1;
			}
			if (k == 0){
				bcb[m] = m - 1;
				fcb[m] = m;
			}
			if (k == HeightNum - 1){
				bcf[m] = m - WidthNum - 1;
				fcf[m] = m - WidthNum;
			}
			/////////////////////////////////////////////////////////
			if ((i == 0) && (k == 0)){
				bcb[m] = m + WidthNum -2;
				fcb[m] = m;
				bcf[m] = m + WidthNum - 2;
			}
			if ((i == WidthNum - 1) && (k == 0)){
				bcb[m] = m - 1;
				fcb[m] = m - WidthNum + 1;
				fcf[m] = m - WidthNum + 1;
			}
			if ((i == 0) && (k == HeightNum - 1)){
				bcb[m] = m - 2;
				bcf[m] = m - 2;
				fcf[m] = m - WidthNum;
			}
			if ((i == WidthNum - 1) && (k == HeightNum - 1)){
				fcb[m] = m - WidthNum - WidthNum + 1;
				bcf[m] = m - WidthNum - 1;
				fcf[m] = m - WidthNum - WidthNum + 1;
			}
			/////////////////////////////////////////////////////////
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition(){

	//�l�p�d�ɂ�z�u����
	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum + i;
			ElePerm[m] = spElePerm * AirRelElePerm;

			//�d�ɏ���ݒ�
			for (int e = 0; e < cylinder_Num; e++){

				double dx = abs(i - cylinder[e].pos_X);
				double dz = abs(k - cylinder[e].pos_Z);
				double distance = sqrt((dx * dx) + (dz * dz));

				//�d�ʂ�ݒ�
				if (distance <= cylinder[e].radius){

					CalcFlag[m] = false;
					if (e < 2){
						Potential[m] = Volt;
					}
					else
					{
						Potential[m] = -Volt;
					}
				}

				//�U�d���ݒ�
				if (distance <= cylinder[e].diele_radius){
					ElePerm[m] = spElePerm * CylinderElePerm;
				}
			}

			//�d�ɏ���ݒ�
			for (int e = 0; e < wallNum_electrode; e++){

				if ((i >= wall_electrode[e].pos_X[0]) && (i <= wall_electrode[e].pos_X[1])
					&& (k >= wall_electrode[e].pos_Z[0]) && (k <= wall_electrode[e].pos_Z[2])){

					CalcFlag[m] = false;
					if (e < 2){
						Potential[m] = Volt;
					}
					else
					{
						Potential[m] = -Volt;
					}				
				}
			}

			//�U�d�̏���ݒ�
			for (int e = 0; e < wallNum_dielectric; e++){

				if ((i >= wall_dielectric[e].pos_X[0]) && (i <= wall_dielectric[e].pos_X[1])
					&& (k >= wall_dielectric[e].pos_Z[0]) && (k <= wall_dielectric[e].pos_Z[2])){
					ElePerm[m] = spElePerm * wallElePerm;
				}
			}
		}
	}
}

void Electric_Calc::Calc_Potential_Laplace(){

	int localCount = 0,m;
	double P, P_fcc, P_bcc, P_ccf, P_ccb;
	double e_bcb, e_fcb, e_bcf, e_fcf;

	#pragma omp single
	cout << "���v���X�̎����v�Z���܂�" << endl;

	do{
		#pragma omp barrier	
		localCount = 0;
				
		#pragma omp single	
		convergentCount = 0;		
		
		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k=k+2){
			for (int i = 0; i<WidthNum; i=i+2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					P_bcc = Potential[bcc[m]];
					P_fcc = Potential[fcc[m]];
					P_ccb = Potential[ccb[m]];
					P_ccf = Potential[ccf[m]];					

					e_bcb = ElePerm[bcb[m]];
					e_fcb = ElePerm[fcb[m]];
					e_bcf = ElePerm[bcf[m]];
					e_fcf = ElePerm[fcf[m]];

					Potential[m] = (1.0 - SOR) * P + SOR * ((P_fcc*(e_fcb + e_fcf) + P_bcc*(e_bcb + e_bcf) + P_ccf*(e_bcf + e_fcf) + P_ccb*(e_bcb + e_fcb)) / (2.0*(e_bcb + e_fcb + e_bcf + e_fcf)));

					if (fabs(Potential[m] - P)>TFC)  localCount++;					
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k = k+2){
			for (int i = 1; i<WidthNum; i = i+2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					P_bcc = Potential[bcc[m]];
					P_fcc = Potential[fcc[m]];
					P_ccb = Potential[ccb[m]];
					P_ccf = Potential[ccf[m]];

					e_bcb = ElePerm[bcb[m]];
					e_fcb = ElePerm[fcb[m]];
					e_bcf = ElePerm[bcf[m]];
					e_fcf = ElePerm[fcf[m]];

					Potential[m] = (1.0 - SOR) * P + SOR * ((P_fcc*(e_fcb + e_fcf) + P_bcc*(e_bcb + e_bcf) + P_ccf*(e_bcf + e_fcf) + P_ccb*(e_bcb + e_fcb)) / (2.0*(e_bcb + e_fcb + e_bcf + e_fcf)));

					if (fabs(Potential[m] - P)>TFC)  localCount++;
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k = k+2){
			for (int i = 1; i<WidthNum; i = i+2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					P_bcc = Potential[bcc[m]];
					P_fcc = Potential[fcc[m]];
					P_ccb = Potential[ccb[m]];
					P_ccf = Potential[ccf[m]];

					e_bcb = ElePerm[bcb[m]];
					e_fcb = ElePerm[fcb[m]];
					e_bcf = ElePerm[bcf[m]];
					e_fcf = ElePerm[fcf[m]];

					Potential[m] = (1.0 - SOR) * P + SOR * ((P_fcc*(e_fcb + e_fcf) + P_bcc*(e_bcb + e_bcf) + P_ccf*(e_bcf + e_fcf) + P_ccb*(e_bcb + e_fcb)) / (2.0*(e_bcb + e_fcb + e_bcf + e_fcf)));

					if (fabs(Potential[m] - P)>TFC)  localCount++;
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k = k+2){
			for (int i = 0; i<WidthNum; i = i+2){

				//�ߓ_�ԍ������߂�
				m = k * WidthNum + i;

				if (CalcFlag[m]){

					//�ߓ_�̓d��	
					P = Potential[m];
					P_bcc = Potential[bcc[m]];
					P_fcc = Potential[fcc[m]];
					P_ccb = Potential[ccb[m]];
					P_ccf = Potential[ccf[m]];

					e_bcb = ElePerm[bcb[m]];
					e_fcb = ElePerm[fcb[m]];
					e_bcf = ElePerm[bcf[m]];
					e_fcf = ElePerm[fcf[m]];

					Potential[m] = (1.0 - SOR) * P + SOR * ((P_fcc*(e_fcb + e_fcf) + P_bcc*(e_bcb + e_bcf) + P_ccf*(e_bcf + e_fcf) + P_ccb*(e_bcb + e_fcb)) / (2.0*(e_bcb + e_fcb + e_bcf + e_fcf)));

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

void Electric_Calc::Calc_E(){

	//�d�E���v�Z����
	if (Type_boundary_X == 0)		Calc_E_X_n();
	else if (Type_boundary_X == 1)	Calc_E_X_r();

	Calc_E_Z_n();
	
	cout << "�d�E���v�Z" << endl;
}

void Electric_Calc::Calc_E_X_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;

			if (i == 0){
				E_x[m] = 0.0;
			}
			else if (i == WidthNum - 1){
				E_x[m] = 0.0;
			}
			else{
				E_x[m] = -(Potential[m + 1] - Potential[m - 1]) / (ElmSize * 2.0);
			}
		}
	}
}

void Electric_Calc::Calc_E_X_r(){

	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;

			if (i == 0){
				E_x[m] = -(Potential[m + 1] - Potential[m + WidthNum - 1]) / (ElmSize * 2.0);
			}
			else if (i == WidthNum - 1){
				E_x[m] = -(Potential[m - WidthNum + 1] - Potential[m - 1]) / (ElmSize * 2.0);
			}
			else{
				E_x[m] = -(Potential[m + 1] - Potential[m - 1]) / (ElmSize * 2.0);
			}
		}
	}
}
void Electric_Calc::Calc_E_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int i = 0; i < WidthNum; i++){

			int m = k * WidthNum + i;

			if (k == 0){
				E_z[m] = 0.0;
			}
			else if (k == HeightNum - 1){
				E_z[m] = 0.0;
			}
			else{
				E_z[m] = -(Potential[m + WidthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
			}
		}
	}
}

void Electric_Calc::OutPut_Potential(char *filename_1, char *filename_2){

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

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			m = k * WidthNum + i;
			fout_potential2 << Potential[m] << '\n';
		}
	}

	fout_potential1.close();
	fout_potential2.close();
}

void Electric_Calc::OutPut_E(char *filename_1, char *filename_2){

	cout << "�d�E�f�[�^���o�͂��܂�" << endl;

	ofstream fout_E1;
	ofstream fout_E2;

	fout_E1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_E2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_E1
		<< "# AVS field file\n" << "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			m = k * WidthNum + i;
			fout_E2 << sqrt(pow(E_x[m], 2) + pow(E_z[m], 2)) << '\n';

		}
	}

	fout_E1.close();
	fout_E2.close();
}

void Electric_Calc::OutPut_EleFieldBinary(char *filename){

	cout << "�����v�Z�p�̓d�E�f�[�^���o�͂��܂�" << endl;

	ofstream ofs;
	ofs.open(filename, ios_base::out | ios_base::binary | ios_base::trunc);
	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){
			m = k * WidthNum + i;
			ofs.write((const char*)&Potential[m], sizeof(double));
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){
			m = k * WidthNum + i;
			ofs.write((const char*)&E_x[m], sizeof(double));
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){
			m = k * WidthNum + i;
			ofs.write((const char*)&E_z[m], sizeof(double));
		}
	}
}

void Electric_Calc::OutPut_Condition(char *filename){

	cout << "�v�Z�������o�͂��܂�" << endl;

	ofstream fout_E;
	fout_E.open(filename, ios_base::out | ios_base::trunc);

	fout_E << "mesh_size =  " << ElmSize << " [m] " << endl;
	fout_E << "x_direction =  " << areaWidth << " [m] " << endl;
	fout_E << "z_direction =  " << areaHeight << " [m] " << endl;
	fout_E << "x_mesh =  " << WidthNum << endl;
	fout_E << "z_mesh =  " << HeightNum << endl;
	fout_E << "Volt =  " << Volt << " [V] " << endl;

	if (Type_boundary_X == 0)		fout_E << "X���� : �m�C�}�����E����" << endl;
	else if (Type_boundary_X == 1)	fout_E << "X���� : �������E����" << endl;

	fout_E << "Z���� : �f�B���N�����E����" << endl;
}
