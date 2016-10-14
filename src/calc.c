#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include"calc.h"
#include"valuables.h"
#include"utilities.h"

//valuables frem runCEA.C
extern T_CEA CEA[3]; // parameter for chamber, throat, exit

//global valuables just in this file
static float At,Dt,Ae,De,Ln,Ac;
//At: throat area, Dt:throat Diam, Ae: exit area, De: exit diameter, Ln: nozzle length, Ac:chamber area
static float Vc,Dc,Lc,Ct;
//Vc: chamber volume, Rc:chamber radius, Lc:chamber length, Ct:chamber thickness

static float mt,mf, mo;
//mt: total fuel mass, mf:fuel mass, mo:mass lox

static float hg,rg,hf,rf,hm,rm,rt;

void calc_chamber_spec(void){

	calc_nozzle();
	calc_chamber();
	calc_chamber_strength();
	calc_fuel_consumption();
	calc_regene();

}

void calc_nozzle(void){

	At = THRUST / (PC * CEA[nozzle_exit].Cf);
	Dt = sqrt(4 * At / M_PI);
	Ae = At * CEA[nozzle_exit].AeAt;
	De = sqrt(4 * Ae / M_PI);

	printf("At = %f[mm^2]\n",At);
	printf("Dt = %f[mm]\n",Dt);
	printf("Ae = %f[mm^2]\n",Ae);
	printf("De = %f[mm]\n",De);

	//calculate conical nozzle length

	Ln = 1.5 * Dt/2 * sin(convert_to(rad,15)) 
			+ (Dt/2 * CEA[nozzle_exit].AeAt - Dt/2)/tan(convert_to(rad,15));

	printf("Ln = %f[mm]\n",Ln);
}

void calc_chamber(void){

	char temp[64];

	Vc = At * convert_to(mm,L_STAR);
	
	printf("Vc = %f[mm^3]\n",Vc);

	while(1){
		printf("input chamber diameter in mm. when satisfied, input type ok\n");
		printf("Dc = ");
		scanf("%s",temp);
		if(strcmp(temp,"ok") == 0){
			break;	
		}else{
			Dc = atof(temp);
			if(Dc == 0){
				printf("error input correct value\n");
			}else{
				Lc = Vc / (Dc * Dc * M_PI / 4);
				Ac = Dc * Dc * M_PI / 4;
				printf("Lc = %f[mm]\n",Lc);
				printf("Ac = %f[mm^2]\n",Ac);
			}
		}
	}
}

void calc_chamber_strength(void){
	Ct = (PC * Dc) / ((2 * TENSILE_STRENGTH / SAFETY_COEFF) - 1.2 * PC);

	printf("Ct = %f\n",Ct);
}

void calc_fuel_consumption(void){
	mt = THRUST / (CEA[nozzle_exit].Isp);
	mf = mt / (OF + 1);
	mo = mt * OF / (OF + 1);
	printf("mt = %f[kg/sec]\n",mt);
	printf("mf = %f[kg/sec]\n",mf);
	printf("mo = %f[kg/sec]\n",mo);
}

void calc_regene(void){

	calc_gas_heat_stansfer_coeff();
	calc_fuel_heat_stansfer_coeff();
	calc_metal_heat_stansfer_coeff();
	calc_total_heat_stansfer_coeff();
	calc_delta_fuel_temp();
}

// using bartz
// needs to be modified for geometrical change in chamber
void calc_gas_heat_stansfer_coeff(void){
	int i = 0;
	float bartz[6];
	hg = 1;

	bartz[0] = 0.026 / pow(convert_to(m,Dt),0.2);
	bartz[1] = pow(CEA[chamber].Visc_gas/10000,0.2) * CEA[chamber].Cp * 1000
				/ pow(CEA[chamber].Prandtl_gas,0.6); 
				//visc changed from mmpoise to kgm/sec, Cp transforms to J
	bartz[2] = pow(convert_to(Pa,PC) / CEA[nozzle_exit].Cstar,0.8);
	bartz[3] = pow(Dt / (1.5 * Dt / 2),0.1);	//if the nozzle are conical
	bartz[4] = pow(At / Ac,0.9); //change here for geometrical parameter.

	
	for(i = 0; i < 5; i++){
		hg *= bartz[i];
//		printf("bartz%d = %f\n",i,bartz[i]);
	}

	rg = 1 / hg;
	printf("hg = %f[W/m^2K]\n",hg);
	printf("rg = %f[m^2K/W]\n",rg);
}


void calc_fuel_heat_stansfer_coeff(void){
	float path_width,path_height,path_area,path_num,hydraulic_diam;
	char temp1[64],temp2[64];

	while(1){
		printf("input fuel path dimention for cooling\n");
		printf("path width [mm] = ");
		scanf("%s",temp1);
		printf("path height[mm] = ");
		scanf("%s",temp2);

		path_width = atof(temp1);
		path_height = atof(temp2);
		
		path_area = path_width * path_height;
		//this hydraulic diameter is only for square type path
		hydraulic_diam = 2 * path_area / ( path_width + path_height);

		printf("path area = %f\n",path_area);
		printf("hydraulic diameter = %f\n",hydraulic_diam);

		printf("if you are satisfied, type ok\n");
		scanf("%s",temp1);

		if(strcmp(temp1,"ok") == 0){
			break;	
		}
	}



	

}

void calc_metal_heat_stansfer_coeff(void){

}

void calc_total_heat_stansfer_coeff(void){

}

void calc_delta_fuel_temp(void){

}

void calc_fuel_cost(void){

}
