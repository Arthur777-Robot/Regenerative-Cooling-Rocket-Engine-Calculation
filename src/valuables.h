#ifndef INCLUDED
#define INCLUDED

/**********************************/
/************Defines***************/
/**********************************/

#define THRUST 5000		//[N] 
#define PC	1.6			//[MPa]
#define OF  2.1			//[ratio]
#define TYPE_FUEL "C10H21,n-decyl"	//write in CEA format
#define TYPE_OXIDIZER "O2(L)"	//write in CEA format
#define DENSITY_FUEL 0.83	//[kg/L] default value is Ethanol
#define DENSITY_OXIDIZER 1.14	//[kg/L] default value is LOX
#define L_STAR 1.4		//[m]
#define TENSILE_STRENGTH 150 //[N/mm^2] 150[N/mm^2] Cu@250deg_c
#define SAFETY_COEFF 2.5
#define CP_F 2093 		//Fuel thermal capacity [J/kgK] 
#define VISC_F 0.00074291		//fuel viscocity at high temp [Kgm/sec]
#define THRM_COND_FUEL 0.133		// fuel thermal conductivity at high temp [W/mK]
#define THRM_COND_METAL 390		// chamber material heat conductivity [W/mK]

typedef struct{
	float Tc;		//Chamber Temp[K]
	float AeAt;		//Throat ratio
	float Cstar;	//[m/sec]
	float Isp;		//[sec]
	float Cf;		//
	float Cp;		//[J/kg*K]
	float Prandtl_gas;	//
	float Visc_gas;	//[kg*m/sec]
	float Mach;
}T_CEA;

typedef enum{
	chamber,
	throat,
	nozzle_exit,
}T_place;

typedef enum{
	mm,
	mm2,
	mm3,
	m,
	m2,
	m3,
	gram,
	kg,
	Pa,
	MPa,
	deg_c,
	deg_k,
	rad,
	deg,
	L,
}T_prefix;


/*****************************/
/*****Physical properties*****/
/*****************************/
#define MOL_L 22.5	//[L/mol]
#define P_ATM 1		//[kg/cm^2]
#define GAS_CONST	//[J/mol*K]


#endif
