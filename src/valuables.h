#ifndef INCLUDED
#define INCLUDED

/**********************************/
/************Defines***************/
/**********************************/

#define THRUST 5000		//[N] 
#define PC	2.0			//[MPa]
#define OF  2.3			//[ratio]
//#define TYPE_FUEL "C10H21,n-decyl"	//write in CEA format
#define TYPE_FUEL "Jet-A(L)"	//write in CEA format
#define TYPE_OXIDIZER "O2(L)"	//write in CEA format
#define DENSITY_FUEL 0.79	//[kg/L] default value is Ethanol
#define DENSITY_OXIDIZER 1.14	//[kg/L] default value is LOX
#define L_STAR 1.4		//[m]
#define TENSILE_STRENGTH 150 //[N/mm^2] 150[N/mm^2] Cu@250deg_c
#define SAFETY_COEFF 2.5
#define CP_F 2093.0 		//Fuel thermal capacity [J/kgK] 
#define VISC_F 0.00074291		//fuel viscocity at high temp [Kgm/sec]
#define THRM_COND_FUEL 0.133		// fuel thermal conductivity at high temp [W/mK]
#define THRM_COND_METAL 390.0		// chamber material heat conductivity [W/mK]

typedef struct{
	double Tc;		//Chamber Temp[K]
	double AeAt;		//Throat ratio
	double Cstar;	//[m/sec]
	double Isp;		//[sec]
	double Cf;		//
	double Cp;		//[J/kg*K]
	double Prandtl_gas;	//
	double Visc_gas;	//[kg*m/sec]
	double Mach;
	double Gamma;
	double Ivac;		//[sec]
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
