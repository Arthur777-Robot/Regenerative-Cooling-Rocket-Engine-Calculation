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
static float Vc,Dc,Lc,Ct,vf;
//Vc: chamber volume, Rc:chamber radius, Lc:chamber length, Ct:chamber thickness

static float mt,mf, mo;
//mt: total fuel mass, mf:fuel mass, mo:mass lox

static float hg,rg,hf,rf,hm,rm,rc,rt;	//rc for carbon heat resist

static float Q, Tcw,Tf_out;

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

	printf("Nozzle throat area: At = %f[mm^2]\n",At);
	printf("Nozzle throat diameter: Dt = %f[mm]\n",Dt);
	printf("Nozzle exit area: Ae = %f[mm^2]\n",Ae);
	printf("Nozzle exit diameter: De = %f[mm]\n",De);

	//calculate conical nozzle length

	Ln = 1.5 * Dt/2 * sin(convert_to(rad,15)) 
			+ (Dt/2 * CEA[nozzle_exit].AeAt - Dt/2)/tan(convert_to(rad,15));

	printf("Nozzle Length from throat to exit: Ln = %f[mm]\n",Ln);
}

void calc_chamber(void){

	char temp[64];

	Vc = At * convert_to(mm,L_STAR);
	
	printf("Chamber volume: Vc = %f[mm^3]\n",Vc);

	while(1){
		printf("input chamber diameter in mm. when satisfied, input type ok\n");
		printf("Chamber Diameter: Dc = ");
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
				printf("Chamber length: Lc = %f[mm]\n",Lc);
				printf("Chamber area: Ac = %f[mm^2]\n",Ac);
			}
		}
	}
}

void calc_chamber_strength(void){
	Ct = (PC * Dc) / ((2 * TENSILE_STRENGTH / SAFETY_COEFF) - 1.2 * PC);

	printf("Chamber Thickness: Ct = %f\n",Ct);
}

void calc_fuel_consumption(void){
	mt = THRUST / (CEA[nozzle_exit].Isp);
	mf = mt / (OF + 1);
	mo = mt * OF / (OF + 1);
	printf("Total Fuel Consumption: mt = %f[kg/sec]\n",mt);
	printf("Fuel Consumption: mf = %f[kg/sec]\n",mf);
	printf("Oxidizer consumption: mo = %f[kg/sec]\n",mo);
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
	printf("Gas heat conductivity: = %f[W/m^2K]\n",hg);
	printf("Gas heat resistance:  = %f[m^2K/W]\n",rg);
}


void calc_fuel_heat_stansfer_coeff(void){ //helical coil type
	float path_width,path_height,path_area,path_num,hydraulic_diam;
	float nyuf,delta_p;
	float Re,Pr,Nu;
	char temp1[64],temp2[64],temp3[64];

	while(1){
		printf("input fuel path dimention for cooling\n");
		printf("Fuel velocity in cooling channel should be from 9-12m/sec\n");
		printf("path width [mm] = ");
		scanf("%s",temp1);
		printf("path height[mm] = ");
		scanf("%s",temp2);
		printf("number of channel = ");
		scanf("%s",temp3);

		path_width = atof(temp1);
		path_height = atof(temp2);
		path_num = atof(temp3);
		
		path_area = path_width * path_height;
		//this hydraulic diameter is only for square type path
		hydraulic_diam = 2 * path_area / ( path_width + path_height);
		vf = mf / DENSITY_FUEL / 1000 / convert_to(m2,path_area) / path_num;

		delta_p = 0.5 * DENSITY_FUEL * 1000 * pow(vf,2);

		printf("Path area = %f[mm^2]\n",path_area);
		printf("Hydraulic diameter = %f[mm]\n",hydraulic_diam);
		printf("Fuel velocity in cooling channel = %f[m/sec]\n",vf);
		printf("Channel pressure drop = %f[MPa]\n",convert_to(MPa,delta_p));

		printf("if you are satisfied, type ok\n");
		scanf("%s",temp1);

		if(strcmp(temp1,"ok") == 0){
			break;	
		}
	}

	nyuf = VISC_F / DENSITY_FUEL / 1000;
	Re = vf * convert_to(m,hydraulic_diam) / nyuf ;
	Pr = VISC_F * CP_F / THRM_COND_FUEL;
	Nu = 0.023 * pow(Re,0.8) * pow(Pr,0.6);

	hf = Nu * THRM_COND_FUEL / convert_to(m,hydraulic_diam);
	rf = 1 / hf; 

	printf("Coefficient of kinematic viscosity = %e\n",nyuf);
	printf("Reynolds number of liquid = %f\n",Re);
	printf("Prandtle number of liquid = %f\n",Pr);
	printf("Nusselt number of liquid = %f\n",Nu);
	printf("Fuel heat conductivity: hf = %f[W/m^2K]\n",hf);
	printf("Fuel heat resistance: rf = %e[m^2K/W]\n",rf);
	
}

void calc_metal_heat_stansfer_coeff(void){

	hm = THRM_COND_METAL;
	rm = convert_to(m,Ct) / hm;

	printf("chamber material heat conductivity = %f\n",hm);
	printf("chamber material heat resistance = %e\n",rm);

}

void calc_total_heat_stansfer_coeff(void){
	float Tcwc;

	rc = 0.000407663;	//carbon heat resistivity

	rt = rg + rf + rm + rc;
	Q = (CEA[chamber].Tc - 298) / rt;
	Tcw = CEA[chamber].Tc - Q * (rg + rc);
	Tcwc = CEA[chamber].Tc - Q * (rg + rc + rm);

	printf("Total heat resistance = %f[m^2K/W]\n",rt);
	printf("Total heat transfer = %f[MW/m^2]\n",Q/1000000);
	printf("Chamber wall temp = %f[deg_c]\n",convert_to(deg_c,Tcw));
	printf("Chamber wall channel temp = %f[deg_c]\n",convert_to(deg_c,Tcwc));
	printf("Check = %f[K]\n",CEA[chamber].Tc - Q * (rt));

}

void calc_delta_fuel_temp(void){
	float A_ct;		//Chamber total area

	A_ct = M_PI * (Dc + 2 * Ct) * Lc;
	A_ct = A_ct * 1.2;		//nozzle are typically 10% of chamber

	Tf_out = convert_to(m2,A_ct) * Q / (vf * CP_F);

	printf("Chamber area total = %f[mm^2]\n",convert_to(m2,A_ct));
	printf("Total heat to wall = %f[W]\n",convert_to(m2,A_ct*Q));
	printf("Fuel delta temp = %f[deg_c]\n",Tf_out);

}

void calc_fuel_cost(void){

}


