#include "Ele_Calc.h"

Electric_Calc::Electric_Calc() :PI(3.14159265358979323846264338328){		//�R���X�g���N�^�D�ϐ��̏����������Ă���D

	//��{�p�����[�^
	ElmSize = 200.0e-6;			//���b�V���T�C�Y
	areaWidth = 50000.0e-6;		//�w�����v�Z�̈�
	areaDepth = 50000.0e-6;		//Y�����̌v�Z�̈�
	areaHeight = 50000.0e-6;	//�y�����v�Z�̈�	
	Type_boundary_X = 0;		//���̋��E�����̐ݒ�@0 = �m�C�}�� / 1 = ���� ����I��
	Type_boundary_Y = 0;		//���̋��E�����̐ݒ�@0 = �m�C�}�� / 1 = ���� ����I��

	spElePerm = 8.854187816e-12;//�U�d��
	AirRelElePerm = 1.00059;	//��Ԃ̔�U�d��
	Volt = 5000;				//����d��		
	TFC = 1.0e-3;				//��������
	SOR = 1.8;					//�����W��	
	OMP_CALC = true;			//���񉻂��邩
	thread = 6;					//����
	unsteady_calc = false;		//����̌v�Z��
	load_ptcl_Poisson = false;	//���q�f�[�^��ǂݍ��ނ��@�|�A�\����
	
	//����d�E�v�Z���s���ꍇ�̂�
	if (unsteady_calc){
		total_T = 0.5;					// �V�~�����[�V��������								
		delta_T = 0.5e-3;				// ���U�����Ԃ͎��萔�� 1/10 �ȉ��ɂȂ�悤��
		// (ParticleRelElePerm * spElePerm)/ conductivity)		-6 : 2.21e-5, -7 : 2.21e-4, -8 : 2.21e-3, -9 : 2.21e-2, -10 : 2.21e-1,	-11 : 2.21,   -12 : 2.21e+1		���萔
		frequency = 1;					//���g��
		phase = 2;						//�d�ʂ̑���
		rise_time = 2.5e-3;				//�d�ʂ̗����オ�莞��
		WriteAVS = 50;					//����̏������݃X�e�b�v
		total_step = (int)(total_T / delta_T);
		phase_No = 2;
		phase_change_flag = true;
		span = 0.0;
	}

	//���q���l������E�|�A�\���̏ꍇ�̂�
	if (load_ptcl_Poisson){
		ParticleRelElePerm = 2.5;		//���q�̔�U�d��
		conductivity = 1.0e-9;			//���d��
	}	
	
	//���d�ɂ�����ꍇ
	cylinder_Xdir_Num = 6 ;					//�d�ɐ�
	cylinder_Xdir.resize(cylinder_Xdir_Num);		
	CylinderElePerm = 5.0;				//�d�ɔ핞�̗U�d��

	for (int i = 0; i < cylinder_Xdir_Num; i++){

		cylinder_Xdir[i].Pos_X = 7500.0e-6;
		cylinder_Xdir[i].Pos_Y = 12500.0e-6 + 5000.0e-6 * i;
		cylinder_Xdir[i].Pos_Z = 9800.0e-6;
		cylinder_Xdir[i].pos_X = (int)(cylinder_Xdir[i].Pos_X / ElmSize + 0.5);
		cylinder_Xdir[i].pos_Y = (int)(cylinder_Xdir[i].Pos_Y / ElmSize + 0.5);
		cylinder_Xdir[i].pos_Z = (int)(cylinder_Xdir[i].Pos_Z / ElmSize + 0.5);

		cylinder_Xdir[i].Radius = 400.0e-6;					//�d�ɔ��a
		cylinder_Xdir[i].Length = 35000.0e-6;				//�d�ɒ���
		cylinder_Xdir[i].Dielectric_Radius = 900.0e-6;		//�d�ɔ핞���܂߂����a
		cylinder_Xdir[i].radius = (int)(cylinder_Xdir[i].Radius / ElmSize + 0.5);
		cylinder_Xdir[i].length = (int)(cylinder_Xdir[i].Length / ElmSize + 0.5);
		cylinder_Xdir[i].diele_radius = (int)(cylinder_Xdir[i].Dielectric_Radius / ElmSize + 0.5);
	}

	//���d�ɂ�����ꍇ
	cylinder_Ydir_Num = 6;					//�d�ɐ�
	cylinder_Ydir.resize(cylinder_Ydir_Num);

	for (int i = 0; i < cylinder_Ydir_Num; i++){

		cylinder_Ydir[i].Pos_X = 12500.0e-6 + 5000.0e-6 * i;
		cylinder_Ydir[i].Pos_Y = 7500.0e-6;
		cylinder_Ydir[i].Pos_Z = 6800.0e-6;
		cylinder_Ydir[i].pos_X = (int)(cylinder_Ydir[i].Pos_X / ElmSize + 0.5);
		cylinder_Ydir[i].pos_Y = (int)(cylinder_Ydir[i].Pos_Y / ElmSize + 0.5);
		cylinder_Ydir[i].pos_Z = (int)(cylinder_Ydir[i].Pos_Z / ElmSize + 0.5);

		cylinder_Ydir[i].Radius = 400.0e-6;				//�d�ɔ��a
		cylinder_Ydir[i].Length = 35000.0e-6;			//�d�ɒ���
		cylinder_Ydir[i].Dielectric_Radius = 900.0e-6;  //�d�ɔ핞���܂߂����a
		cylinder_Ydir[i].radius = (int)(cylinder_Ydir[i].Radius / ElmSize + 0.5);
		cylinder_Ydir[i].length = (int)(cylinder_Ydir[i].Length / ElmSize + 0.5);
		cylinder_Ydir[i].diele_radius = (int)(cylinder_Ydir[i].Dielectric_Radius / ElmSize + 0.5);
	}	

	//�~��������ꍇ
	hollow_cylNum = 1;	
	hollow_cylinder.resize(hollow_cylNum);
	hollow_cylinderElePerm = 3.1;

	for (int i = 0; i < hollow_cylNum; i++){

		hollow_cylinder[i].Pos_X = 25000.0e-6;
		hollow_cylinder[i].Pos_Y = 25000.0e-6;
		hollow_cylinder[i].Pos_Z = 0.0e-6;
		hollow_cylinder[i].pos_X = (int)(hollow_cylinder[i].Pos_X / ElmSize + 0.5);
		hollow_cylinder[i].pos_Y = (int)(hollow_cylinder[i].Pos_Y / ElmSize + 0.5);
		hollow_cylinder[i].pos_Z = (int)(hollow_cylinder[i].Pos_Z / ElmSize + 0.5);

		hollow_cylinder[i].Inside_Radius = 17500.0e-6;
		hollow_cylinder[i].Outside_Radius = 25000.0e-6;
		hollow_cylinder[i].Length = 16000.0e-6;				
		hollow_cylinder[i].inside_radius = (int)(hollow_cylinder[i].Inside_Radius / ElmSize + 0.5);
		hollow_cylinder[i].outside_radius = (int)(hollow_cylinder[i].Outside_Radius / ElmSize + 0.5);
		hollow_cylinder[i].length = (int)(hollow_cylinder[i].Length / ElmSize + 0.5);
	}

	//���������ꍇ
	wallNum = 1;
	wall.resize(wallNum);
	wallElePerm = 2.5;

	for (int i = 0; i < wallNum; i++){

		wall[i].Pos_X[0] = 0;
		wall[i].Pos_Y[0] = 0;
		wall[i].Pos_Z[0] = 0;

		wall[i].Pos_X[1] = areaWidth;
		wall[i].Pos_Y[1] = 0;
		wall[i].Pos_Z[1] = 0;

		wall[i].Pos_X[2] = 0;
		wall[i].Pos_Y[2] = 0;
		wall[i].Pos_Z[2] = 5900.0e-6;

		wall[i].Pos_X[3] = areaWidth;
		wall[i].Pos_Y[3] = 0;
		wall[i].Pos_Z[3] = 5900.0e-6;		

		for (int n = 0; n < 4; n++){
			wall[i].pos_X[n] = (int)(wall[i].Pos_X[n] / ElmSize + 0.5);
			wall[i].pos_Y[n] = (int)(wall[i].Pos_Y[n] / ElmSize + 0.5);
			wall[i].pos_Z[n] = (int)(wall[i].Pos_Z[n] / ElmSize + 0.5);
		}

		wall[i].Length = areaDepth;
		wall[i].length = (int)(wall[i].Length / ElmSize + 0.5);		
	}

	//�ȍ~�͕ϐ��̏������E�v�Z�p�����[�^�̗��U��
	//��{�p�����[�^�Ɋւ���
	WidthNum = (int)(areaWidth / ElmSize + 0.5);
	DepthNum = (int)(areaDepth / ElmSize + 0.5);
	HeightNum = (int)(areaHeight / ElmSize + 0.5);
	cellNum_2D = WidthNum * HeightNum;
	cellNum_3D = WidthNum * DepthNum * HeightNum;
	if (cellNum_3D == 0)	cellNum_3D = 1;
	convergentCount = 0;

	ElePerm.assign(cellNum_3D, 0.0);	Potential.assign(cellNum_3D, 0.0);
	E_x.assign(cellNum_3D, 0.0);	E_y.assign(cellNum_3D, 0.0);	E_z.assign(cellNum_3D, 0.0);
	bcc.assign(cellNum_3D, 0);		fcc.assign(cellNum_3D, 0);		cbc.assign(cellNum_3D, 0);
	cfc.assign(cellNum_3D, 0);		ccb.assign(cellNum_3D, 0);		ccf.assign(cellNum_3D, 0);
	bbb.assign(cellNum_3D, 0);		bbf.assign(cellNum_3D, 0);		bfb.assign(cellNum_3D, 0);		bff.assign(cellNum_3D, 0);
	fbb.assign(cellNum_3D, 0);		fbf.assign(cellNum_3D, 0);		ffb.assign(cellNum_3D, 0);		fff.assign(cellNum_3D, 0);

	CalcFlag = new bool[cellNum_3D];
	for (int i = 0; i < cellNum_3D; i++)	CalcFlag[i] = true;

	//���q���l������ꍇ�E�|�A�\���̎��̏ꍇ
	if (load_ptcl_Poisson){
		ptclNum = 0;
		ptcl_num_node.assign(cellNum_3D, -1);
		ptcl_num_ele.assign(cellNum_3D, -1);
		Conductivity.assign(cellNum_3D, 0);
		Charge.assign(cellNum_3D, 0);
	}

	cout << "�v�Z������ݒ芮��" << endl;
}

Electric_Calc::~Electric_Calc(){

	if (CalcFlag !=0)	delete[] CalcFlag;
}

void Electric_Calc::Setting_Boundary_Condition_V(){

	//���̓d�ʂɊւ��鋫�E�����ݒ�
	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				bcc[m] = m - 1;
				fcc[m] = m + 1;
				cbc[m] = m - WidthNum;
				cfc[m] = m + WidthNum;
				ccb[m] = m - (WidthNum * DepthNum);
				ccf[m] = m + (WidthNum * DepthNum);
			}
		}
	}

	if (Type_boundary_X == 0)			Setting_Boundary_Condition_V_X_n();
	else if (Type_boundary_X == 1)		Setting_Boundary_Condition_V_X_l();
	
	if (Type_boundary_Y == 0)			Setting_Boundary_Condition_V_Y_n();
	else if (Type_boundary_Y == 1)		Setting_Boundary_Condition_V_Y_l();
	
	Setting_Boundary_Condition_V_Z_n();
	
	cout << "�d�ʂɊւ���[�̋��E�����ݒ�" << endl;
}

void Electric_Calc::Setting_Boundary_Condition_V_X_n(){

	for(int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (i == 0)					bcc[m] = m + 1;
				if (i == WidthNum - 1)		fcc[m] = m - 1;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_X_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (i == 0)					bcc[m] = m + WidthNum - 1;
				if (i == WidthNum - 1)		fcc[m] = m - WidthNum + 1;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_Y_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (j == 0)					cbc[m] = m + WidthNum;
				if (j == DepthNum - 1)		cfc[m] = m - WidthNum;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_Y_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (j == 0)					cbc[m] = m + (HeightNum * WidthNum) - WidthNum;
				if (j == DepthNum - 1)		cfc[m] = m - (HeightNum * WidthNum) + WidthNum;
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_V_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				if (k == 0)					ccb[m] = m + (WidthNum * DepthNum);
				if (k == HeightNum - 1)		ccf[m] = m - (WidthNum * DepthNum);
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e(){

	//���̗U�d���Ɋւ��鋫�E�����ݒ�
	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				bbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
				fbb[m] = m - (WidthNum*DepthNum) - WidthNum;
				bfb[m] = m - (WidthNum*DepthNum) - 1;
				ffb[m] = m - (WidthNum*DepthNum);
				bbf[m] = m - WidthNum - 1;
				fbf[m] = m - WidthNum;
				bff[m] = m - 1;
				fff[m] = m;
			}
		}
	}
	Setting_Boundary_Condition_e_X_n_corner();          //�ȉ��[�Ɗp�̗U�d���Ɋւ��Ă͌v�Z�ɑ傫�ȉe����^���Ȃ����߁A�m�C�}���Ƃ��Ĉ���(�d�ʂɊւ��Ă͍l�����ׂ�)
	Setting_Boundary_Condition_e_Y_n_corner();          //�ȉ��[�Ɗp�̗U�d���Ɋւ��Ă͌v�Z�ɑ傫�ȉe����^���Ȃ����߁A�m�C�}���Ƃ��Ĉ���
	Setting_Boundary_Condition_e_Z_n_corner();          //�ȉ��[�Ɗp�̗U�d���Ɋւ��Ă͌v�Z�ɑ傫�ȉe����^���Ȃ����߁A�m�C�}���Ƃ��Ĉ���
	Setting_Boundary_Condition_e_X_n_Y_n_corner();      //�ȉ��[�Ɗp�̗U�d���Ɋւ��Ă͌v�Z�ɑ傫�ȉe����^���Ȃ����߁A�m�C�}���Ƃ��Ĉ���
	Setting_Boundary_Condition_e_X_n_Z_n_corner();      //�ȉ��[�Ɗp�̗U�d���Ɋւ��Ă͌v�Z�ɑ傫�ȉe����^���Ȃ����߁A�m�C�}���Ƃ��Ĉ���
	Setting_Boundary_Condition_e_Y_n_Z_n_corner();      //�ȉ��[�Ɗp�̗U�d���Ɋւ��Ă͌v�Z�ɑ傫�ȉe����^���Ȃ����߁A�m�C�}���Ƃ��Ĉ���
	Setting_Boundary_Condition_e_X_n_Y_n_Z_n_corner(); 

	cout << "�U�d���Ɋւ���[�̋��E�����ݒ�" << endl;
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (i == 0){
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bfb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - WidthNum;
					bff[m] = m;
				}
				if (i == WidthNum - 1){
					fbb[m] = m - (WidthNum * DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum * DepthNum) - 1;
					fbf[m] = m - WidthNum - 1;
					fff[m] = m - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_Y_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (j == 0){
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - 1;
					fbf[m] = m;
				}
				if (j == DepthNum - 1){
					bfb[m] = m - (WidthNum * DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum * DepthNum) - WidthNum;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (k == 0){
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum;
					bfb[m] = m - 1;
					ffb[m] = m;
				}
				if (k == HeightNum - 1){
					bbf[m] = m - (WidthNum * DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum * DepthNum) - WidthNum;
					bff[m] = m - (WidthNum * DepthNum) - 1;
					fff[m] = m - (WidthNum * DepthNum);
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_Y_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((i == 0) && (j == 0)){
					bbb[m] = m - (WidthNum*DepthNum);
					fbb[m] = m - (WidthNum*DepthNum);
					bfb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m;
					fbf[m] = m;
					bff[m] = m;
				}
				if ((i == 0) && (j == DepthNum - 1)){
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bbf[m] = m - WidthNum;
					bff[m] = m - WidthNum;
					fff[m] = m - WidthNum;
				}
				if ((i == WidthNum - 1) && (j == 0)){
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum) - 1;
					ffb[m] = m - (WidthNum*DepthNum) - 1;
					bbf[m] = m - 1;
					fbf[m] = m - 1;
					fff[m] = m - 1;
				}
				if ((i == WidthNum - 1) && (j == DepthNum - 1)){
					fbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - WidthNum - 1;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((i == 0) && (k == 0)){
					bbb[m] = m - WidthNum;
					fbb[m] = m - WidthNum;
					bfb[m] = m;
					ffb[m] = m;
					bbf[m] = m - WidthNum;					
					bff[m] = m;
				}
				if ((i == 0) && (k == HeightNum - 1)){
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;					
					bfb[m] = m - (WidthNum*DepthNum);					
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					bff[m] = m - (WidthNum*DepthNum);
					fff[m] = m - (WidthNum*DepthNum);
				}
				if ((i == WidthNum - 1) && (k == 0)){
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum - 1;
					bfb[m] = m - 1;
					ffb[m] = m - 1;					
					fbf[m] = m - WidthNum - 1;					
					fff[m] = m - 1;
				}
				if ((i == WidthNum - 1) && (k == HeightNum - 1)){					
					fbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;					
					ffb[m] = m - (WidthNum*DepthNum) - 1;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bff[m] = m - (WidthNum*DepthNum) - 1;
					fff[m] = m - (WidthNum*DepthNum) - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_Y_n_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((j == 0) && (k == 0)){
					bbb[m] = m - 1;
					fbb[m] = m;
					bfb[m] = m - 1;
					ffb[m] = m;
					bbf[m] = m - 1;
					fbf[m] = m;
				}
				if ((j == 0) && (k == HeightNum - 1)){
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - (WidthNum*DepthNum) - 1;
					fbf[m] = m - (WidthNum*DepthNum);
					bff[m] = m - (WidthNum*DepthNum) - 1;
					fff[m] = m - (WidthNum*DepthNum);
				}
				if ((j == DepthNum - 1) && (k == 0)){
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum;
					bfb[m] = m - WidthNum - 1;
					ffb[m] = m - WidthNum;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum;
				}
				if ((j == DepthNum - 1) && (k == HeightNum - 1)){
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					bff[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fff[m] = m - (WidthNum*DepthNum) - WidthNum;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_e_X_n_Y_n_Z_n_corner(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if ((i == 0) && (j == 0) && (k == 0)){//fff
					bbb[m] = m;
					fbb[m] = m;
					bfb[m] = m;
					ffb[m] = m;
					bbf[m] = m;
					fbf[m] = m;
					bff[m] = m;
				}
				if ((i == 0) && (j == 0) && (k == HeightNum - 1)){//ffb
					bbb[m] = m - (WidthNum*DepthNum);
					fbb[m] = m - (WidthNum*DepthNum);
					bfb[m] = m - (WidthNum*DepthNum);
					bbf[m] = m - (WidthNum*DepthNum);
					fbf[m] = m - (WidthNum*DepthNum);
					bff[m] = m - (WidthNum*DepthNum);
					fff[m] = m - (WidthNum*DepthNum);
				}
				if ((i == 0) && (j == DepthNum -1) && (k == 0)){//fbf
					bbb[m] = m - WidthNum;
					fbb[m] = m - WidthNum;
					bfb[m] = m - WidthNum;
					ffb[m] = m - WidthNum;
					bbf[m] = m - WidthNum;
					bff[m] = m - WidthNum;
					fff[m] = m - WidthNum;
				}
				if ((i == 0) && (j == DepthNum - 1) && (k == HeightNum -1)){//fbb
					bbb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum;
					bff[m] = m - (WidthNum*DepthNum) - WidthNum;
					fff[m] = m - (WidthNum*DepthNum) - WidthNum;
				}
				if ((i == WidthNum - 1) && (j == 0) && (k == 0)){//bff
					bbb[m] = m - 1;
					fbb[m] = m - 1;
					bfb[m] = m - 1;
					ffb[m] = m - 1;
					bbf[m] = m - 1;
					fbf[m] = m - 1;
					fff[m] = m - 1;
				}
				if ((i == WidthNum - 1) && (j == 0) && (k == HeightNum - 1)){//bfb
					bbb[m] = m - (WidthNum*DepthNum) - 1;
					fbb[m] = m - (WidthNum*DepthNum) - 1;
					ffb[m] = m - (WidthNum*DepthNum) - 1;
					bbf[m] = m - (WidthNum*DepthNum) - 1;
					fbf[m] = m - (WidthNum*DepthNum) - 1;
					bff[m] = m - (WidthNum*DepthNum) - 1;
					fff[m] = m - (WidthNum*DepthNum) - 1;
				}
				if ((i == WidthNum - 1) && (j == DepthNum - 1) && (k == 0)){//bbf
					bbb[m] = m - WidthNum - 1;
					fbb[m] = m - WidthNum - 1;
					bfb[m] = m - WidthNum - 1;
					ffb[m] = m - WidthNum - 1;
					fbf[m] = m - WidthNum - 1;
					bff[m] = m - WidthNum - 1;
					fff[m] = m - WidthNum - 1;
				}
				if ((i == WidthNum - 1) && (j == DepthNum - 1) && (k == HeightNum - 1)){//bbb
					fbb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bfb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					ffb[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fbf[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					bff[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
					fff[m] = m - (WidthNum*DepthNum) - WidthNum - 1;
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_1(){

	cout << "�d�ɓ��̋��E�������Z�b�e�C���O" << endl;

	//���E�����d�ʃZ�b�g
	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				ElePerm[m] = spElePerm * AirRelElePerm;			

				//�d�ɏ���ݒ�
				for (int e = 0; e < cylinder_Xdir_Num; e++){

					double dy = abs(j - cylinder_Xdir[e].pos_Y);
					double dz = abs(k - cylinder_Xdir[e].pos_Z);
					double distance = sqrt((dy * dy) + (dz * dz));

					//�d�ʂ�ݒ�
					if ((distance <= cylinder_Xdir[e].radius) && (cylinder_Xdir[e].pos_X <= i) && (i <= (cylinder_Xdir[e].pos_X + cylinder_Xdir[e].length))){
						cylinder_Xdir[e].cyl_mesh.push_back(m);
						CalcFlag[m] = false;
						Potential[m] = - Volt;
					}

					//�U�d���ݒ�
					if ((distance <= cylinder_Xdir[e].diele_radius) && (cylinder_Xdir[e].pos_X <= i) && (i <= (cylinder_Xdir[e].pos_X + cylinder_Xdir[e].length))){
						ElePerm[m] = spElePerm * CylinderElePerm;
					}
				}

				//�d�ɏ���ݒ�
				for (int e = 0; e < cylinder_Ydir_Num; e++){

					double dx = abs(i - cylinder_Ydir[e].pos_X);
					double dz = abs(k - cylinder_Ydir[e].pos_Z);
					double distance = sqrt((dx * dx) + (dz * dz));

					//�d�ʂ�ݒ�
					if ((distance <= cylinder_Ydir[e].radius) && (cylinder_Ydir[e].pos_Y <= j) && (j <= (cylinder_Ydir[e].pos_Y + cylinder_Ydir[e].length))){
						cylinder_Ydir[e].cyl_mesh.push_back(m);
						CalcFlag[m] = false;
						Potential[m] = Volt;
					}

					//�U�d���ݒ�
					if ((distance <= cylinder_Ydir[e].diele_radius) && (cylinder_Ydir[e].pos_Y <= j) && (j <= (cylinder_Ydir[e].pos_Y + cylinder_Ydir[e].length))){
						ElePerm[m] = spElePerm * CylinderElePerm;
					}
				}

				//����~������ݒ�
				for (int e = 0; e < hollow_cylNum; e++){

					double dx = abs(i - hollow_cylinder[e].pos_X);
					double dy = abs(j - hollow_cylinder[e].pos_Y);
					double distance = sqrt((dx * dx) + (dy * dy));

					//�U�d���ݒ�
					if ((hollow_cylinder[e].inside_radius <= distance)
						&& (distance <= hollow_cylinder[e].outside_radius) 
						&& (hollow_cylinder[e].pos_Z <= k) 
						&& (k <= (hollow_cylinder[e].pos_Z + hollow_cylinder[e].length))){

						hollow_cylinder[e].hollow_cyl_mesh.push_back(m);
						ElePerm[m] = spElePerm * hollow_cylinderElePerm;
					}
				}

				//�������ݒ�
				for (int e = 0; e < wallNum; e++){

					//�U�d���ݒ�
					if ((i >= wall[e].pos_X[0]) && (i <= wall[e].pos_X[1])
						&& (j >= wall[e].pos_Y[0]) && (j <= wall[e].pos_Y[0] + wall[e].length)
						&& (k >= wall[e].pos_Z[0]) && (k <= wall[e].pos_Z[2])){

						wall[e].wall_mesh.push_back(m);
						ElePerm[m] = spElePerm * wallElePerm;
					}
				}
			}
		}
	}
}

void Electric_Calc::Setting_Boundary_Condition_2(int step, char *filename){

	//���E�����d�ʃZ�b�g
	//�d�ʗ����オ�莞
	if (((step * delta_T) >= span) && ((step * delta_T) <= span + rise_time)){

		//�d�E��؂�ւ��邩
		if (phase_change_flag){

			if (phase_No == 1)			phase_No = 2;
			else if (phase_No == 2)		phase_No = 1;			
			phase_change_flag = false;
		}

		if (phase_No == 1){			
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = -Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}		
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}
		}
		else if (phase_No == 2){			
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = Potential[m] = Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = -Volt*(1 - exp(-((step * delta_T) - span) / (rise_time / 4.0)));
				}
			}
		}

		//�d�ʗ����オ�莞�Ԃ𒴂����u��
	}
	else if ((step * delta_T) > span + rise_time){

		phase_change_flag = true;
		span += (1.0 / frequency) / phase;		

		if (phase_No == 1){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = -Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = Volt;
				}
			}
		}
		else if (phase_No == 2){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = -Volt;
				}
			}
		}

		//�d�ʗ����オ�莞�Ԃ𒴂�����
	}
	else{
		if (phase_No == 1){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = -Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = Volt;
				}
			}
		}
		else if (phase_No == 2){
			for (int e = 0; e < cylinder_Xdir_Num; e++){
				for (int m = 0; m < cylinder_Xdir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Xdir[e].cyl_mesh[m]] = Volt;
				}
			}
			for (int e = 0; e < cylinder_Ydir_Num; e++){
				for (int m = 0; m < cylinder_Ydir[e].cyl_mesh.size(); m++){
					Potential[cylinder_Ydir[e].cyl_mesh[m]] = -Volt;
				}
			}
		}
	}

	//�d�ʃf�[�^���o��	
	ofstream fout;
	if (step == 0){
		fout.open(filename, ios_base::out, ios_base::trunc);
	}
	else{
		fout.open(filename, ios_base::app);
	}

	for (int e = 0; e < cylinder_Xdir_Num; e++){
		fout << Potential[cylinder_Xdir[e].cyl_mesh[0]] << ",";
	}

	for (int e = 0; e < cylinder_Ydir_Num; e++){
		fout << Potential[cylinder_Ydir[e].cyl_mesh[0]] << ",";
	}
	fout << endl;
	fout.close();
}

void Electric_Calc::Ptcl_Data(char *filename){

	//���q�f�[�^�̃I�[�v��
	cout << "���q�f�[�^��ǂݍ���ł��܂�..." << endl;
	ifstream fin(filename);

	//�G���[�`�F�b�N
	if (!fin){
		cout << "Error.Can't open 1ptcl data file." << endl;
		exit(1);
	}

	string temp1, temp2;
	int p = 0;

	//���q�f�[�^�̓ǂ݂���
	//1�s�ڂ̐��l��ǂݍ��݁A���q����ݒ�
	getline(fin, temp1);
	p = int(temp1.find(","));
	temp2 = temp1.substr(0, p);
	ptclNum = int_string(temp2);

	P_C.resize(ptclNum);

	for (int i = 0; i < ptclNum; i++){
		P_C[i].Ptcl_Pos_X = 0.0;
		P_C[i].Ptcl_Pos_Y = 0.0;
		P_C[i].Ptcl_Pos_Z = 0.0;
		P_C[i].Ptcl_Radius = 0.0;
		P_C[i].Ptcl_Charge = 0.0;
		P_C[i].sp_G = 0.0;
	}	
	
	//2�s�ڂ�ǂݍ���ŁA��������
	getline(fin, temp1);

	//�e���q�̃p�����[�^��ݒ�
	for (int i = 0; i < ptclNum; i++){

		//��s�܂�܂�ǂݍ���
		getline(fin, temp1);

		//split�����s�@getline�œǂݍ��񂾕������","�ŋ�؂�Astr_List�ɕۑ�
		list<string> str_List = split(temp1, ",");

		//�C�e���[�^���擾  ","�ŋ�؂�ꂽ�f�[�^�̒��ŁA1�ԍŏ��̂��̂ɃC�e���[�^�����킹��
		list<string>::iterator iter = str_List.begin();

		//�C�e���[�^���w������f�[�^��������
		temp2 = *iter;
		P_C[i].Ptcl_Pos_X = double_string(temp2);
		++iter;

		//�C�e���[�^��i�߂āA���̃f�[�^��������@�ȍ~�͌J��Ԃ�
		temp2 = *iter;
		P_C[i].Ptcl_Pos_Y = double_string(temp2);
		++iter;

		temp2 = *iter;
		P_C[i].Ptcl_Pos_Z = double_string(temp2);
		++iter;

		++iter;
		++iter;
		++iter;
		++iter;
		++iter;
		++iter;

		temp2 = *iter;
		P_C[i].sp_G = double_string(temp2);
		++iter;

		temp2 = *iter;
		P_C[i].Ptcl_Radius = 0.5 * double_string(temp2);
		++iter;

		temp2 = *iter;
		P_C[i].Ptcl_Charge = 4.0 / 3.0 * PI * pow(P_C[i].Ptcl_Radius, 3) * P_C[i].sp_G * double_string(temp2) * 1.0E-3;
	}

	fin.close();
}

void Electric_Calc::Particle_Node_Mapping(){

	cout << "�ߓ_��ɗ��q��ݒu���܂�" << endl;

	int P_x, P_y, P_z, P_r, m;
	double rel_x, rel_y, rel_z, rel;

	//���q���S�ʒu�����
	for (int p = 0; p < ptclNum; p++){

		//���q�ʒu�𐮐����W��
		P_x = (int)(P_C[p].Ptcl_Pos_X / ElmSize + 0.5);
		P_y = (int)(P_C[p].Ptcl_Pos_Y / ElmSize + 0.5);
		P_z = (int)(P_C[p].Ptcl_Pos_Z / ElmSize + 0.5);
		P_r = (int)(P_C[p].Ptcl_Radius / ElmSize + 0.5);

		for (int k = 0; k <= HeightNum; k++){
			for (int j = P_y - P_r - 1; j <= P_y + P_r + 1; j++){
				for (int i = P_x - P_r - 1; i <= P_x + P_r + 1; i++){

					//���q�ƑΏۃZ���̑��΋������Z�o
					rel_x = (P_x - i) * ElmSize;
					rel_y = (P_y - j) * ElmSize;
					rel_z = (P_z - k) * ElmSize;
					rel = rel_x*rel_x + rel_y*rel_y + rel_z*rel_z;

					//�d�Ȃ��Ă���Η��q�����Z�b�g
					if (rel == 0){
						m = k * WidthNum * DepthNum + j * WidthNum + i;
						ptcl_num_node[m] = p;
						P_C[p].node.push_back(m);
					}
					else{
						if (rel <= pow(P_C[p].Ptcl_Radius, 2)){

							//�ߓ_�̔ԍ������߂�
							//X������-�̐ߓ_
							if (i<0){
								//Y������-�̐ߓ_
								if (j<0){
									m = k * WidthNum * DepthNum + j * WidthNum + i + WidthNum * DepthNum + WidthNum;
									//Y�������ő�l�𒴂����ߓ_
								}
								else if (j >= DepthNum){
									m = k * WidthNum * DepthNum + j * WidthNum + i - WidthNum * DepthNum + WidthNum;
								}
								else{
									m = k * WidthNum * DepthNum + j * WidthNum + i + WidthNum;
								}
								//X�������ő�l�𒴂����ߓ_
							}
							else if (i >= WidthNum){
								//Y������-�̐ߓ_
								if (j<0){
									m = k * WidthNum * DepthNum + j * WidthNum + i + WidthNum * DepthNum - WidthNum;
									//Y�������ő�l�𒴂����ߓ_
								}
								else if (j >= DepthNum){
									m = k * WidthNum * DepthNum + j * WidthNum + i - WidthNum * DepthNum - WidthNum;
								}
								else{
									m = k * WidthNum * DepthNum + j * WidthNum + i - WidthNum;
								}
							}
							else{
								//Y������-�̐ߓ_
								if (j<0){
									m = k * WidthNum * DepthNum + j * WidthNum + WidthNum * DepthNum + i;
									//Y�������ő�l�𒴂����ߓ_
								}
								else if (j >= DepthNum){
									m = k * WidthNum * DepthNum + j * WidthNum - WidthNum * DepthNum + i;
								}
								else{
									m = k * WidthNum * DepthNum + j * WidthNum + i;
								}
							}

							ptcl_num_node[m] = p;
							P_C[p].node.push_back(m);
						}
					}
				}
			}
		}
	}

	/*
	
	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open("check.fld", ios_base::out | ios_base::trunc);
	fout_charge2.open("check.dat", ios_base::out | ios_base::trunc);

	fout_charge1
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << "check.dat" << ' ' << "filetype=ascii\n";
	
	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_charge2 << ptcl_num_node[m] << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();

	*/

}

void Electric_Calc::Particle_Ele_Mapping(){

	cout << "�v�f��ɗ��q��ݒu���܂�" << endl;

	int p_node_0, p_node_1, p_node_2, p_node_3, p_node_4, p_node_5, p_node_6, p_node_7;

	//���q���S�ʒu�����
	for (int p = 0; p < ptclNum; p++){
		int node = (int)P_C[p].node.size();

		for (int i = 0; i < node; i++){
			int m = P_C[p].node[i];

			p_node_0 = ptcl_num_node[m];
			p_node_1 = ptcl_num_node[m + 1];
			p_node_2 = ptcl_num_node[m + WidthNum];
			p_node_3 = ptcl_num_node[m + WidthNum + 1];
			p_node_4 = ptcl_num_node[m + WidthNum*DepthNum];
			p_node_5 = ptcl_num_node[m + WidthNum*DepthNum + 1];
			p_node_6 = ptcl_num_node[m + WidthNum*DepthNum + WidthNum];
			p_node_7 = ptcl_num_node[m + WidthNum*DepthNum + WidthNum + 1];

			//�Ώۂ̗v�f���\������S�ߓ_�ɁA�������q�̏�񂪋L������Ă���΁A���q�Ɨv�f���d�Ȃ��Ă���
			if ((p_node_0 == p_node_1) && (p_node_0 == p_node_2) && (p_node_0 == p_node_3) && (p_node_0 == p_node_4) && (p_node_0 == p_node_5) && (p_node_0 == p_node_6) && (p_node_0 == p_node_7)){

				//�v�f�̔ԍ����擾���A���q���Ɨ��q�̗U�d�����L������
				ElePerm[m] = spElePerm * ParticleRelElePerm;
				Conductivity[m] = conductivity;
				ptcl_num_ele[m] = ptcl_num_node[m];
				P_C[p].ele.push_back(m);
			}
		}
	}

	/*
	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open("check.fld", ios_base::out | ios_base::trunc);
	fout_charge2.open("check.dat", ios_base::out | ios_base::trunc);

	fout_charge1
	<< "# AVS field file\n" << "ndim = 3\n"
	<< "dim1 = " << WidthNum << '\n'
	<< "dim2 = " << DepthNum << '\n'
	<< "dim3 = " << HeightNum << '\n'
	<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
	<< "field = uniform\n" << '\n'
	<< "variable 1 file=" << "check.dat" << ' ' << "filetype=ascii\n";

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){
				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_charge2 << ptcl_num_ele[m] << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();
	*/

}

void Electric_Calc::Particle_count_charge(){

	cout << "���q�̍ŊO�ʒu�ɑ��݂���ߓ_����肵�܂�" << endl;

	int n_0, e_1, e_2, e_3, e_4, e_5, e_6;

	//���q���\������O���ߓ_�ɂ̂ݓd�ׂ�z�u����
	for (int p = 0; p < ptclNum; p++){
		int ele = (int)P_C[p].ele.size();

		for (int d = 0; d < ele; d++){
			int m = P_C[p].ele[d];
			int i = (m % (WidthNum*DepthNum)) % WidthNum;
			int j = (m % (WidthNum*DepthNum)) / WidthNum;
			int k = m / (WidthNum*DepthNum);

			n_0 = ptcl_num_ele[m];
			e_5 = m + WidthNum*DepthNum;
			e_6 = m - WidthNum*DepthNum;

			//���ڗv�f��X�EY�����ɗאڂ���v�f�����擾
			//X������0�̗v�f
			if (i == 0){
				//Y������0�̗v�f
				if (j == 0){
					e_1 = m + 1;
					e_2 = m + WidthNum - 1;
					e_3 = m + WidthNum;
					e_4 = m + WidthNum*DepthNum - WidthNum;
					//Y������MAX�̗v�f
				}
				else if (j == DepthNum - 1){
					e_1 = m + 1;
					e_2 = m + WidthNum - 1;
					e_3 = m - WidthNum*DepthNum + WidthNum;
					e_4 = m - WidthNum;
					//Y�������[�v�f�ȊO
				}
				else{
					e_1 = m + 1;
					e_2 = m + WidthNum - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum;
				}
				//X������MAX�̗v�f
			}
			else if (i == WidthNum - 1){
				//Y������0�̗v�f
				if (j == 0){
					e_1 = m + 1 - WidthNum;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum + WidthNum*DepthNum;
					//Y������MAX�̗v�f
				}
				else if (j == DepthNum - 1){
					e_1 = m - WidthNum + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum - WidthNum*DepthNum;
					e_4 = m - WidthNum;
					//Y�������[�v�f�ȊO
				}
				else{
					e_1 = m + 1 - WidthNum;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum;
				}
				//X�������[�v�f�ȊO
			}
			else{
				//Y������0�̗v�f
				if (j == 0){
					e_1 = m + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m + WidthNum*DepthNum - WidthNum;
					//Y������MAX
				}
				else if (j == DepthNum - 1){
					e_1 = m + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum - WidthNum*DepthNum;
					e_4 = m - WidthNum;
					//Y�������[�_�ȊO
				}
				else{
					e_1 = m + 1;
					e_2 = m - 1;
					e_3 = m + WidthNum;
					e_4 = m - WidthNum;
				}
			}

			//���ڗv�f�Ǝ��͗v�f�ɓ������q�ԍ����ݒ肳��Ă���Ηv�f�͗אڂ��Ă���B����āA�ŊO�ʒu�ł͂Ȃ� �ȉ��͍ŊO�ʒu�����߂Ă���
			if (!(n_0 == ptcl_num_ele[e_1])){
				P_C[n_0].contain.push_back(m + 1);
				P_C[n_0].contain.push_back(m + 1 + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1 + WidthNum);
			}
			if (!(n_0 == ptcl_num_ele[e_2])){
				P_C[n_0].contain.push_back(m);
				P_C[n_0].contain.push_back(m + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + WidthNum);
			}
			if (!(n_0 == ptcl_num_ele[e_3])){
				P_C[n_0].contain.push_back(m + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum + 1);
				P_C[n_0].contain.push_back(m + WidthNum + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum + WidthNum*DepthNum + 1);
			}
			if (!(n_0 == ptcl_num_ele[e_4])){
				P_C[n_0].contain.push_back(m);
				P_C[n_0].contain.push_back(m + 1);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1);
			}
			if (!(n_0 == ptcl_num_ele[e_5])){
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + 1);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum*DepthNum + WidthNum + 1);
			}

			//Z�����̉��[�ʂ̗v�f
			if (k == 0){

				P_C[n_0].contain.push_back(m);
				P_C[n_0].contain.push_back(m + 1);
				P_C[n_0].contain.push_back(m + WidthNum);
				P_C[n_0].contain.push_back(m + WidthNum + 1);

				//Z�����̒[�ʈȊO�̗v�f
			}
			else{
				if (!(n_0 == ptcl_num_ele[e_6])){
					P_C[n_0].contain.push_back(m);
					P_C[n_0].contain.push_back(m + 1);
					P_C[n_0].contain.push_back(m + WidthNum);
					P_C[n_0].contain.push_back(m + WidthNum + 1);
				}
			}
		}
	}

	//�����v�f���폜����
	for (int i = 0; i<ptclNum; i++){

		//sort���ėv�f�̏��Ԃ������ŕ��ёւ���
		sort(P_C[i].contain.begin(), P_C[i].contain.end());
		//unique�����ĘA�����铯���v�f��A�������A���c�����S�~��erase�ō폜����
		P_C[i].contain.erase(unique(P_C[i].contain.begin(), P_C[i].contain.end()), P_C[i].contain.end());
	}

	/*
	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open("check.fld", ios_base::out | ios_base::trunc);
	fout_charge2.open("check.dat", ios_base::out | ios_base::trunc);

	fout_charge1
	<< "# AVS field file\n" << "ndim = 3\n"
	<< "dim1 = " << WidthNum << '\n'
	<< "dim2 = " << DepthNum << '\n'
	<< "dim3 = " << HeightNum << '\n'
	<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
	<< "field = uniform\n" << '\n'
	<< "variable 1 file=" << "check.dat" << ' ' << "filetype=ascii\n";

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				int temp = 0;

				//sort���ėv�f�̏��Ԃ������ŕ��ёւ���
				for (int n = 0; n < P_C[0].contain.size(); n++)		if (P_C[0].contain[n] == m)		temp = 10;
				fout_charge2 << temp << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();
	*/

}

void Electric_Calc::Particle_charge_Mapping(){

	//���q�̓d�ׂ�z�u����
	for (int i = 0; i<ptclNum; i++){

		int Size = int(P_C[i].contain.size());

		for (int j = 0; j<Size; j++){

			int cell_num = P_C[i].contain[j];
			Charge[cell_num] = P_C[i].Ptcl_Charge / Size;
		}
	}
}

void Electric_Calc::Particle_clear(){

	ptcl_num_node.clear();
	ptcl_num_ele.clear();
}


void Electric_Calc::Calc_Potential_Laplace(){

	int localCount = 0;
	double P, P_fcc, P_bcc, P_ccf, P_ccb, P_cfc, P_cbc;
	double e_bbb, e_bbf, e_bfb, e_bff, e_fbb, e_fbf, e_ffb, e_fff;
	int m;

	do{
		#pragma omp barrier
		localCount = 0;

		#pragma omp single
		convergentCount = 0;

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k+=2){
			for (int j = 0; j<DepthNum; j+=2){
				for (int i = 0; i<WidthNum; i+=2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];							
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR 
										* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
											+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
											+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)) 
												/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));
											
						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}					
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

		#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
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

void Electric_Calc::Calc_Potential_Poisson(){

	int localCount = 0;
	double Q, P, P_fcc, P_bcc, P_ccf, P_ccb, P_cfc, P_cbc;
	double e_bbb, e_bbf, e_bfb, e_bff, e_fbb, e_fbf, e_ffb, e_fff;
	int m;

	do{
#pragma omp barrier
		localCount = 0;

#pragma omp single
		convergentCount = 0;

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 0; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for nowait
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 0; j<DepthNum; j += 2){
				for (int i = 0; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp for
		for (int k = 1; k<HeightNum; k += 2){
			for (int j = 1; j<DepthNum; j += 2){
				for (int i = 1; i<WidthNum; i += 2){

					//�ߓ_�ԍ������߂�
					m = k * WidthNum * DepthNum + j * WidthNum + i;

					if (CalcFlag[m]){

						//�ߓ_�̓d��	
						Q = Charge[m];
						P = Potential[m];
						P_bcc = Potential[bcc[m]];
						P_fcc = Potential[fcc[m]];
						P_cbc = Potential[cbc[m]];
						P_cfc = Potential[cfc[m]];
						P_ccb = Potential[ccb[m]];
						P_ccf = Potential[ccf[m]];

						e_bbb = ElePerm[bbb[m]];
						e_fbb = ElePerm[fbb[m]];
						e_bfb = ElePerm[bfb[m]];
						e_ffb = ElePerm[ffb[m]];
						e_bbf = ElePerm[bbf[m]];
						e_fbf = ElePerm[fbf[m]];
						e_bff = ElePerm[bff[m]];
						e_fff = ElePerm[fff[m]];

						Potential[m] = (1.0 - SOR) * P + SOR
							* ((P_bcc*(e_bbb + e_bfb + e_bbf + e_bff) + P_fcc*(e_fbb + e_ffb + e_fbf + e_fff)
							+ P_cbc*(e_bbb + e_fbb + e_bbf + e_fbf) + P_cfc*(e_bfb + e_ffb + e_bff + e_fff)
							+ P_ccb*(e_bbb + e_fbb + e_bfb + e_ffb) + P_ccf*(e_bbf + e_fbf + e_bff + e_fff)
							+ (4 * Q / ElmSize))
							/ (3.0*(e_bbb + e_fbb + e_bfb + e_ffb + e_bbf + e_fbf + e_bff + e_fff)));

						if (fabs(Potential[m] - P)>TFC)  localCount++;
					}
				}
			}
		}

#pragma omp critical
		convergentCount += localCount;

#pragma omp barrier

	} while (convergentCount);
}

void Electric_Calc::Calc_E(){

	//�d�E���v�Z����
	if (Type_boundary_X == 0)		Calc_E_X_n();
	else if (Type_boundary_X == 1)	Calc_E_X_l();

	if (Type_boundary_Y == 0)		Calc_E_Y_n();
	else if (Type_boundary_Y == 1)	Calc_E_Y_l();

	Calc_E_Z_n();
}

void Electric_Calc::Calc_E_X_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

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
}

void Electric_Calc::Calc_E_X_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (i == 0){
					E_x[m] = -(Potential[m + 1] - Potential[m + WidthNum - 2]) / (ElmSize * 2.0);
				}
				else if (i == WidthNum - 1){
					E_x[m] = -(Potential[m - WidthNum + 2] - Potential[m - 1]) / (ElmSize * 2.0);
				}
				else{
					E_x[m] = -(Potential[m + 1] - Potential[m - 1]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Calc_E_Y_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (j == 0){
					E_y[m] = 0.0;
				}
				else if (j == DepthNum - 1){
					E_y[m] = 0.0;
				}
				else{
					E_y[m] = -(Potential[m + WidthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Calc_E_Y_l(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (j == 0){
					E_y[m] = -(Potential[m + WidthNum] - Potential[m + (WidthNum - 2)*DepthNum]) / (ElmSize * 2.0);
				}
				else if (j == DepthNum - 1){
					E_y[m] = -(Potential[m - (WidthNum - 2)*DepthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
				}
				else{
					E_y[m] = -(Potential[m + WidthNum] - Potential[m - WidthNum]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Calc_E_Z_n(){

	for (int k = 0; k < HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){

				int m = k * WidthNum * DepthNum + j * WidthNum + i;

				if (k == 0){
					E_z[m] = 0.0;
				}
				else if (k == HeightNum - 1){
					E_z[m] = 0.0;
				}
				else{
					E_z[m] = -(Potential[m + (WidthNum*DepthNum)] - Potential[m - (WidthNum*DepthNum)]) / (ElmSize * 2.0);
				}
			}
		}
	}
}

void Electric_Calc::Ini_P(){

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				//�ߓ_�ԍ������߂�
				int m = k * WidthNum * DepthNum + j * WidthNum + i;
				Potential[m] = 0.0;
			}
		}
	}
}

void Electric_Calc::Calc_Ohm(){

	int localCount = 0;
	int n_1, n_2, n_3, n_4, n_5, n_6;
	double e_1, e_2, e_3, e_4, e_5, e_6, e_7, e_8;
	double Jxfwd, Jxbck, Jyfwd, Jybck, Jzfwd, Jzbck;		//�e�����̓d�����x

	//�d�ׂ̌v�Z
	for (int p = 0; p < ptclNum; p++){
		int node = (int)P_C[p].node.size();

		for (int i = 0; i < node; i++){

			int m = P_C[p].node[i];
							
				n_1 = m + 1;
				n_2 = m - 1;
				n_3 = m + WidthNum;
				n_4 = m - WidthNum;
				n_5 = m + WidthNum*DepthNum;
				n_6 = m - WidthNum*DepthNum;

				e_1 = Conductivity[m - WidthNum*DepthNum - WidthNum - 1];
				e_2 = Conductivity[m - WidthNum*DepthNum - WidthNum];
				e_3 = Conductivity[m - WidthNum*DepthNum - 1];
				e_4 = Conductivity[m - WidthNum*DepthNum];
				e_5 = Conductivity[m - WidthNum - 1];
				e_6 = Conductivity[m - WidthNum];
				e_7 = Conductivity[m - 1];
				e_8 = Conductivity[m];

				Jxfwd = -1.0 * (e_2 + e_4 + e_6 + e_8) / 4 * (Potential[n_1] - Potential[m]) / ElmSize;
				Jxbck = -1.0 * (e_1 + e_3 + e_5 + e_7) / 4 * (Potential[m] - Potential[n_2]) / ElmSize;
				Jyfwd = -1.0 * (e_3 + e_4 + e_7 + e_8) / 4 * (Potential[n_3] - Potential[m]) / ElmSize;
				Jybck = -1.0 * (e_1 + e_2 + e_5 + e_6) / 4 * (Potential[m] - Potential[n_4]) / ElmSize;
				Jzfwd = -1.0 * (e_5 + e_6 + e_7 + e_8) / 4 * (Potential[n_5] - Potential[m]) / ElmSize;
				Jzbck = -1.0 * (e_1 + e_2 + e_3 + e_4) / 4 * (Potential[m] - Potential[n_6]) / ElmSize;

				Charge[m] += (-Jxfwd + Jxbck - Jyfwd + Jybck - Jzfwd + Jzbck) * ElmSize * ElmSize * delta_T;
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
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_potential2 << Potential[m] << '\n';
			}
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
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){

				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_E2 << sqrt(pow(E_x[m], 2) + pow(E_y[m], 2) + pow(E_z[m], 2)) << '\n';
			}
		}
	}

	fout_E1.close();
	fout_E2.close();
}

void Electric_Calc::OutPut_Charge(char *filename_1, char *filename_2){

	cout << "�d�ו��z�f�[�^���o�͂��܂�" << endl;

	ofstream fout_charge1;
	ofstream fout_charge2;

	fout_charge1.open(filename_1, ios_base::out | ios_base::trunc);
	fout_charge2.open(filename_2, ios_base::out | ios_base::trunc);

	fout_charge1
		<< "# AVS field file\n" << "ndim = 3\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << DepthNum << '\n'
		<< "dim3 = " << HeightNum << '\n'
		<< "nspace = 3\n" << "veclen = 1\n" << "data = double\n"
		<< "field = uniform\n" << '\n'
		<< "variable 1 file=" << filename_2 << ' ' << "filetype=ascii\n";

	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i<WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				fout_charge2 << Charge[m] * 10e+12 << '\n';
			}
		}
	}

	fout_charge1.close();
	fout_charge2.close();
}

void Electric_Calc::OutPut_Ptcl(char *filename){

	cout << "���q�f�[�^���o�͂��܂�" << endl;
	ofstream fout;
	fout.open(filename, ios_base::out | ios_base::trunc);
	fout << "# Micro AVS Geom:2.00\n";
	fout << "sphere" << endl;
	fout << "Movement" << endl;
	fout << "color" << endl;
	fout << ptclNum << endl;

	for (int i = 0; i < ptclNum; i++){
		fout << P_C[i].Ptcl_Pos_X << " " << P_C[i].Ptcl_Pos_Y << " " << P_C[i].Ptcl_Pos_Z << " " << P_C[i].Ptcl_Radius << " 1 0 0" << endl;
	}

	fout << "disjoint polygon" << endl;
	fout << "floor" << endl;
	fout << "facet" << endl;
	fout << "color" << endl;
	fout << 1 << endl;
	double gray_color = 0.8;

	fout << "4" << endl;
	fout << "0" << " " << "0" << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
	fout << areaWidth << " " << "0" << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
	fout << areaWidth << " " << areaDepth << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
	fout << "0" << " " << areaDepth << " " << "0" << " " << gray_color << " " << gray_color << " " << gray_color << " " << endl;
}

void Electric_Calc::OutPut_Force_1(char *filename){

	ofstream fout;
	fout.open(filename, ios_base::out | ios_base::trunc);

	int P_x_1 = (int)(P_C[0].Ptcl_Pos_X / ElmSize + 0.5);
	int P_y_1 = (int)(P_C[0].Ptcl_Pos_Y / ElmSize + 0.5);
	int P_z_1 = (int)(P_C[0].Ptcl_Pos_Z / ElmSize + 0.5);
	int m_1 = P_z_1 * WidthNum * DepthNum + P_y_1 * WidthNum + P_x_1;

	double E_x_1 = E_x[m_1 + 1] * E_x[m_1 + 1] + E_y[m_1 + 1] * E_y[m_1 + 1] + E_z[m_1 + 1] * E_z[m_1 + 1];
	double E_xx_1 = E_x[m_1 - 1] * E_x[m_1 - 1] + E_y[m_1 - 1] * E_y[m_1 - 1] + E_z[m_1 - 1] * E_z[m_1 - 1];
	double E_y_1 = E_x[m_1 + WidthNum] * E_x[m_1 + WidthNum] + E_y[m_1 + WidthNum] * E_y[m_1 + WidthNum] + E_z[m_1 + WidthNum] * E_z[m_1 + WidthNum];
	double E_yy_1 = E_x[m_1 - WidthNum] * E_x[m_1 - WidthNum] + E_y[m_1 - WidthNum] * E_y[m_1 - WidthNum] + E_z[m_1 - WidthNum] * E_z[m_1 - WidthNum];
	double E_z_1 = E_x[m_1 + WidthNum * DepthNum] * E_x[m_1 + WidthNum * DepthNum] + E_y[m_1 + WidthNum * DepthNum] * E_y[m_1 + WidthNum * DepthNum] + E_z[m_1 + WidthNum * DepthNum] * E_z[m_1 + WidthNum * DepthNum];
	double E_zz_1 = E_x[m_1 - WidthNum * DepthNum] * E_x[m_1 - WidthNum * DepthNum] + E_y[m_1 - WidthNum * DepthNum] * E_y[m_1 - WidthNum * DepthNum] + E_z[m_1 - WidthNum * DepthNum] * E_z[m_1 - WidthNum * DepthNum];

	double temp = 2 * PI * spElePerm * AirRelElePerm * (ParticleRelElePerm - AirRelElePerm) / (ParticleRelElePerm - 2 * AirRelElePerm);

	fout << "���q1" << endl;
	fout << "����" << "," << "�O��_X" << "," << "�O��_Y" << "," << "�O��_Z" << endl;
	fout << "0" << "," << P_C[0].Ptcl_Charge * E_x[m_1] << "," << P_C[0].Ptcl_Charge * E_y[m_1] << "," << P_C[0].Ptcl_Charge * E_z[m_1] << endl;

	fout << "���q1" << endl;
	fout << "����" << "," << "�O��_X" << "," << "�O��_Y" << "," << "�O��_Z" << endl;
	fout << "0" << "," << temp * (E_x_1 - E_xx_1) / (2.0*ElmSize) * pow(P_C[0].Ptcl_Radius, 3) << "," << temp * (E_y_1 - E_yy_1) / (2.0*ElmSize) * pow(P_C[0].Ptcl_Radius, 3) << "," << temp * (E_z_1 - E_zz_1) / (2.0*ElmSize) * pow(P_C[0].Ptcl_Radius, 3) << endl;

	fout << "�ȉ��C���[�v" << endl;
	fout << "���q1" << endl;
	fout << "����" << "," << "�O��_X" << "," << "�O��_Y" << "," << "�O��_Z" << endl;

	fout.close();
}

void Electric_Calc::OutPut_Force_2(char *filename, int step){

	ofstream fout;
	fout.open(filename, ios_base::app);

	double Force_1_x = 0, Force_1_y = 0, Force_1_z = 0;

	for (int i = 0; i < P_C[0].node.size(); i++){

		int m = P_C[0].node[i];
		Force_1_x += Charge[m] * E_x[m];
		Force_1_y += Charge[m] * E_y[m];
		Force_1_z += Charge[m] * E_z[m];
	}

	fout << step*delta_T << "," << Force_1_x << "," << Force_1_y << "," << Force_1_z << endl;

	fout.close();
}

void Electric_Calc::OutPut_EleFieldBinary(char *filename){

	cout << "�����v�Z�p�̓d�E�f�[�^���o�͂��܂�" << endl;

	ofstream ofs;
	ofs.open(filename, ios_base::out | ios_base::binary | ios_base::trunc);
	int m;

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&Potential[m], sizeof(double));
			}
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&E_x[m], sizeof(double));
			}
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j < DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&E_y[m], sizeof(double));
			}
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int j = 0; j<DepthNum; j++){
			for (int i = 0; i < WidthNum; i++){
				m = k * WidthNum * DepthNum + j * WidthNum + i;
				ofs.write((const char*)&E_z[m], sizeof(double));
			}
		}
	}
}

void Electric_Calc::GenerateMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat){

	ofstream movie_1;
	ofstream movie_2;
	ofstream movie_3;
	ofstream movie_4;
	ofstream movie_5;
	ofstream movie_6;

	movie_1.open(filename_P_fld, ios_base::out, ios_base::trunc);
	movie_2.open(filename_E_fld, ios_base::out, ios_base::trunc);
	movie_3.open(filename_C_fld, ios_base::out, ios_base::trunc);
	movie_4.open(filename_P_dat, ios_base::out, ios_base::trunc);
	movie_5.open(filename_E_dat, ios_base::out, ios_base::trunc);
	movie_6.open(filename_C_dat, ios_base::out, ios_base::trunc);

	movie_1 << "# AVS field file\n" << "#\n"
		<< "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n"
		<< "veclen = 1\n"
		<< "data = double\n"
		<< "field = uniform\n" << '\n';

	movie_2 << "# AVS field file\n" << "#\n"
		<< "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n"
		<< "veclen = 1\n"
		<< "data = double\n"
		<< "field = uniform\n" << '\n';

	movie_3 << "# AVS field file\n" << "#\n"
		<< "ndim = 2\n"
		<< "dim1 = " << WidthNum << '\n'
		<< "dim2 = " << HeightNum << '\n'
		<< "nspace = 2\n"
		<< "veclen = 1\n"
		<< "data = double\n"
		<< "field = uniform\n" << '\n';

	movie_1.close();
	movie_2.close();
	movie_3.close();
	movie_4.close();
	movie_5.close();
	movie_6.close();
}

void Electric_Calc::WriteMovieFile(char *filename_P_fld, char *filename_P_dat, char *filename_E_fld, char *filename_E_dat, char *filename_C_fld, char *filename_C_dat, int T_step, int Count){

	ofstream movie_1;
	ofstream movie_2;
	ofstream movie_3;
	ofstream movie_4;
	ofstream movie_5;
	ofstream movie_6;

	movie_1.open(filename_P_fld, ios_base::app);
	movie_2.open(filename_E_fld, ios_base::app);
	movie_3.open(filename_C_fld, ios_base::app);
	movie_4.open(filename_P_dat, ios_base::app);
	movie_5.open(filename_E_dat, ios_base::app);
	movie_6.open(filename_C_dat, ios_base::app);

	int STEP;
	STEP = Count*(cellNum_2D + 1);

	int P_y = (int)(P_C[0].Ptcl_Pos_Y / ElmSize + 0.5);

	movie_1 << "time file=" << filename_P_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP << ' ' << "close=1\n"
		<< "variable 1 file=" << filename_P_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP + 1 << ' ' << "close=1\n" << "EOT\n";

	movie_2 << "time file=" << filename_E_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP << ' ' << "close=1\n"
		<< "variable 1 file=" << filename_E_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP + 1 << ' ' << "close=1\n" << "EOT\n";

	movie_3 << "time file=" << filename_C_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP << ' ' << "close=1\n"
		<< "variable 1 file=" << filename_C_dat << ' ' << "filetype=ascii" << ' ' << "skip=" << STEP + 1 << ' ' << "close=1\n" << "EOT\n";

	movie_4 << "time value=" << T_step * delta_T << '\n';
	movie_5 << "time value=" << T_step * delta_T << '\n';
	movie_6 << "time value=" << T_step * delta_T << '\n';

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum * DepthNum + P_y * WidthNum + i;
			movie_4 << Potential[m] << '\n';
		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum * DepthNum + P_y * WidthNum + i;
			movie_5 << sqrt(pow(E_x[m], 2) + pow(E_y[m], 2) + pow(E_z[m], 2)) << '\n';

		}
	}

	for (int k = 0; k<HeightNum; k++){
		for (int i = 0; i<WidthNum; i++){

			int m = k * WidthNum * DepthNum + P_y * WidthNum + i;
			movie_6 << Charge[m] * 10e+12 << '\n';

		}
	}

	movie_1.close();
	movie_2.close();
	movie_3.close();
	movie_4.close();
	movie_5.close();
	movie_6.close();
}

void Electric_Calc::OutPut_Condition(char *filename){

	cout << "�v�Z�������o�͂��܂�" << endl;

	ofstream fout_E;
	fout_E.open(filename, ios_base::out | ios_base::trunc);

	fout_E << "mesh_size =  " << ElmSize << " [m] " << endl;
	fout_E << "x_direction =  " << areaWidth << " [m] " << endl;
	fout_E << "y_direction =  " << areaDepth << " [m] " << endl;
	fout_E << "z_direction =  " << areaHeight << " [m] " << endl;
	fout_E << "x_mesh =  " << WidthNum << endl;
	fout_E << "y_mesh =  " << DepthNum << endl;
	fout_E << "z_mesh =  " << HeightNum << endl;
	fout_E << "Volt =  " << Volt << " [V] " << endl;

	fout_E.close();
}