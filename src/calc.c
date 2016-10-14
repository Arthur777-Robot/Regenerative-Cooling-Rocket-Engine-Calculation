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
void calc_gas_heat_stansfer_coeff(void){
	float bartz1,bartz2,bartz3,bartz4,bartz5;

	bartz1 = 0.026 / pow(convert_to(m,Dt),0.2);
	bartz2 = pow(CEA[chamber].Visc_gas,0.2) * CEA[chamber].Cp * 1000
				/ pow(CEA[chamber].Prandtl_gas,0.6); // Cp transforms to J
	printf("bartz1 = %f\n",bartz1);
	printf("bartz2 = %f\n",bartz2);
}


void calc_fuel_heat_stansfer_coeff(void){

}

void calc_metal_heat_stansfer_coeff(void){

}

void calc_total_heat_stansfer_coeff(void){

}

void calc_delta_fuel_temp(void){

}

void calc_fuel_cost(void){

}
