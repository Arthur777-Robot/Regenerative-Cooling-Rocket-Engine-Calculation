#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include"calc.h"
#include"valuables.h"
#include"utilities.h"
#include"dxf.h"

#define char_num 2000

//valuables frem runCEA.C
extern T_CEA CEA[3]; // parameter for chamber, throat, exit

//global valuables just in this file
static double At,Dt,Ae,De,Ln,Ac;
//At: throat area, Dt:throat Diam, Ae: exit area, De: exit diameter, Ln: nozzle length, Ac:chamber area
static double Vc,Dc,Lc,Ct,vf;
//Vc: chamber volume, Rc:chamber radius, Lc:chamber length, Ct:chamber thickness

static double mt,mf, mo;
//mt: total fuel mass, mf:fuel mass, mo:mass lox

static double hg,rg,hf,rf,hm,rm,rc,rt;	//rc for carbon heat resist

static double Q, Tcw,Tf_out;

static double throat_axis;
static int exit_axis;

static double Chamber_x[char_num];		//chamber geometry
static double Chamber_y[char_num];	//chamber geometry


void calc_chamber_spec(void){

	int i = 0;
	FILE *fp,*gp;
	
	double gas_heat_conductivity[char_num],
		total_heat[char_num],temp_chamber[char_num],temp_channel[char_num],
		fuel_temp[char_num];
	double tmp_Cp,tmp_visc_gas,tmp_prandtl,tmp_gas_temp;

	calc_nozzle();
	calc_chamber();
	calc_chamber_strength();
	calc_fuel_consumption();
	
	calc_fuel_heat_transfer();
	calc_metal_heat_transfer();

	for( i = 0; i < exit_axis;i++){
		tmp_visc_gas = CEA_coeff(i,CEA[chamber].Visc_gas,
				CEA[throat].Visc_gas,CEA[nozzle_exit].Visc_gas);
		tmp_Cp = CEA_coeff(i,CEA[chamber].Cp,CEA[throat].Cp,
				CEA[nozzle_exit].Cp);
		tmp_prandtl = CEA_coeff(i,CEA[chamber].Prandtl_gas,
				CEA[throat].Prandtl_gas,CEA[nozzle_exit].Prandtl_gas);
		tmp_gas_temp = CEA_coeff(i,CEA[chamber].Tc,CEA[throat].Tc,
				CEA[nozzle_exit].Tc);

		calc_gas_heat_transfer(Chamber_y[i],tmp_visc_gas,tmp_Cp,tmp_prandtl);
		calc_total_heat_transfer(i,tmp_gas_temp,total_heat,temp_chamber,
				temp_channel);
	}
	
	fuel_temp[exit_axis + 1] = 0;
	for(i = exit_axis - 1; i >=0; i--){
		fuel_temp[i] = fuel_temp[i + 1] + calc_delta_fuel_temp(Chamber_y[i],total_heat[i]);
	}
	printf("delta temp = %f\n",fuel_temp[10]);

	fp = fopen("data.txt","w");
	for(i = 0; i < exit_axis; i++){
		fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\n",
				Chamber_x[i],Chamber_y[i],temp_chamber[i],
				temp_channel[i],total_heat[i]/1000000,fuel_temp[i]);
	}
	fflush(fp);
	fclose(fp);

	gp = popen("gnuplot -persist","w");
	fprintf(gp,"set datafile separator \"\t\"\n");
	fprintf(gp,"set size square\n");
	fprintf(gp,"set title 'Jet-A 5kN 2Mpa Connical nozzle'\n");
	fprintf(gp,"set key outside\n");
	fprintf(gp,"set ytics nomirror\n");
	fprintf(gp,"set y2tics\n");
	fprintf(gp,"set grid\n");
	fprintf(gp,"set xlabel '[mm]'\n");
	fprintf(gp,"set ylabel '[mm],[deg_c]'\n");
	fprintf(gp,"set y2label '[MW/m^2],[deg_c]'\n");
	fprintf(gp,"set xrange[-10:%f]\n",Chamber_x[exit_axis]+100);
	fprintf(gp,"set yrange[-10:%f]\n",Chamber_x[exit_axis]+100);
	fprintf(gp,"set y2range[0:%d]\n",20);
	fprintf(gp,"plot \"data.txt\" u 1:2 w l t \"chamber geom\"\n");
	fprintf(gp,"replot \"data.txt\" u 1:3 w l t \"chamber wall temp\"\n");
	fprintf(gp,"replot \"data.txt\" u 1:4 w l t \"channel wall temp\"\n");
	fprintf(gp,"replot \"data.txt\" u 1:6 w l t \"Fuel delta temp\"\n");
	fprintf(gp,"replot \"data.txt\" u 1:5 w l t \"Heat Transfer 2nd axis\" axes x1y2\n");
//	fprintf(gp,"replot \"data.txt\" u 1:6 w l t \"Fuel delta temp\"\n");


	pclose(gp);
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
			+ (De/2 - Dt/2)/tan(convert_to(rad,15));

	printf("Nozzle Length from throat to exit: Ln = %f[mm]\n",Ln);
}

void calc_chamber(void){

	char temp[64];

	Vc = At * convert_to(mm,L_STAR);
	
	printf("Chamber volume: Vc = %f[mm^3]\n",Vc);

	while(1){
		printf("input chamber diameter in mm. when satisfied, input type ok\n");
		printf("chamber Diameter should be around %d ~ %d mm\n",(int)Dt*2,(int)Dt*3);
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
				calc_chamber_geom();
				// select nozzletype
				calc_conical_nozzle();
//				calc_foelsch_nozzle();
				plot_chamber();

				//draw chamber geometry in gnuplot
				printf("Chamber length: Lc = %f[mm]\n",Lc);
				printf("Chamber area: Ac = %f[mm^2]\n",Ac);
			}
		}

	}
}

void calc_chamber_strength(void){
	Ct = 1;
//	Ct = (PC * Dc) / ((2 * TENSILE_STRENGTH / SAFETY_COEFF) - 1.2 * PC);

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



// using bartz
// needs to be modified for geometrical change in chamber
void calc_gas_heat_transfer(double chamber_y,double tmp_visc_gas, double tmp_Cp, double tmp_prandtl_gas){
	int j = 0;
	double bartz[6];
	hg = 1;

	bartz[0] = 0.026 / pow(convert_to(m,Dt),0.2);
	bartz[1] = pow(tmp_visc_gas/10000,0.2) * tmp_Cp * 1000
				/ pow(tmp_prandtl_gas,0.6); 
				//visc changed from mmpoise to kgm/sec, Cp transforms to J
	bartz[2] = pow(convert_to(Pa,PC) / CEA[nozzle_exit].Cstar,0.8);
	bartz[3] = pow(Dt / (1.5 * Dt / 2),0.1);	//if the nozzle are conical
	bartz[4] = pow(At / (chamber_y*chamber_y*M_PI),0.9); //change here for geometrical parameter.

	for(j = 0; j < 5; j++){
		hg *= bartz[j];
//		printf("bartz%d = %f\n",j,bartz[j]);
	}

	rg = 1 / hg;
//	printf("Gas heat conductivity: = %f[W/m^2K]\n",hg);
//	printf("Gas heat resistance:  = %f[m^2K/W]\n",rg);
	
}


void calc_fuel_heat_transfer(void){ //helical coil type
	double path_width,path_height,path_area,path_num,hydraulic_diam;
	double nyuf,delta_p;
	double Re,Pr,Nu;
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

void calc_metal_heat_transfer(void){

	hm = THRM_COND_METAL;
	rm = convert_to(m,Ct) / hm;

	printf("chamber material heat conductivity = %f\n",hm);
	printf("chamber material heat resistance = %e\n",rm);

}

void calc_total_heat_transfer(int i,double chamber_temp,double *total_heat,double *temp_chamber,double *temp_channel){
	double Tcwc;

//	rc = 0.000407663;	//carbon heat resistivity
	rc = 0.0;	//carbon heat resistivity

	rt = rg + rf + rm + rc;
	Q = (chamber_temp - 298) / rt;
	Tcw = chamber_temp - Q * (rg + rc);
	Tcwc = chamber_temp - Q * (rg + rc + rm);

	total_heat[i] = Q;
	temp_chamber[i] = convert_to(deg_c,Tcw);
	temp_channel[i] = convert_to(deg_c,Tcwc);

//	printf("Total heat resistance = %f[m^2K/W]\n",rt);
//	printf("Total heat transfer = %f[MW/m^2]\n",total_heat[i]);
//	printf("Chamber wall temp = %f[deg_c]\n",temp_chamber[i]);
//	printf("Chamber wall channel temp = %f[deg_c]\n",temp_channel[i]);
//	printf("rg = %f\trf = %f\trm = %f\trc = %f\n",rg,rf,rm,rc);
//	printf("Check = %f[K]\n",chamber_temp - Q * (rt));

}

double calc_delta_fuel_temp(double chamber_y, double total_heat){
	double A_ct;		//Chamber total area
	double temp;

	A_ct = 2 * M_PI * (chamber_y + Ct);

//	Tf_out = convert_to(m2,A_ct) * total_heat / (vf * CP_F);
	Tf_out = total_heat*convert_to(m2,A_ct)/CP_F/mf;

//	printf("Chamber area total = %f[mm^2]\n",convert_to(m2,A_ct));
//	printf("Total heat to wall = %f[W]\n",convert_to(m2,A_ct*total_heat));
//	printf("Fuel delta temp = %f[deg_c]\n",Tf_out);
//	printf("total_heat = %f\n",total_heat);
//	printf("A_ct = %f\n",A_ct);
//	printf("total heat = %f\n",total_heat);
//	printf("temp = %f\n",Tf_out);
	return Tf_out;
}

void calc_fuel_cost(void){

}

int calc_chamber_geom(void){

	FILE *fp;
	int i = 0;
	int angle1 = 20,size1;
	double x[5],y[5];
	double area[5],total_area;


	size1 = (int)Dt*2;

	//calculate coordinates of chamber
	x[1] = 0;
	y[1] = Dc/2;
	x[2] = size1 * cos(convert_to(rad,90 - angle1)) + x[1];
	y[2] = size1 * sin(convert_to(rad,90 - angle1)) + y[1] - size1;
	y[3] = Dt / 2 * 1.5 * sin(convert_to(rad,270 - angle1)) + Dt / 2 + Dt / 2 * 1.5;
	x[3] = (y[2]-y[3])/tan(convert_to(rad,angle1)) + x[2];
	x[4] = -Dt / 2 * 1.5 * cos(convert_to(rad,270 - angle1)) + x[3];
	y[4] = Dt / 2;
	

	//calculate chamber area by section
	area[1] = pow(size1,2)*M_PI*angle1/360 + (x[2]*(y[2]-y[1]-size1)/2) + (y[1]-size1)*x[2];
	area[2] = (y[2]-y[3])*(x[3]-x[2])/2 + y[3]*(x[3]-x[2]);
	area[3] = (x[4]-x[3])*y[3] - (pow(Dt/2*1.5,2)*M_PI*angle1/360 - (x[4]-x[3])*y[3]/2);

	total_area = area[1]+area[2]+area[3];
	area[0] = Dc*Lc/2 - total_area;

	x[0] = 0; y[0] = Dc/2;

	for(i = 1; i < 5; i++){
		x[i] += area[0]/(Dc/2);
//		printf("x[%d] = %f,y[%d] = %f\n",i,x[i],i,y[i]);
	}

	throat_axis = x[4];

	for(i = 0; i<x[1];i++){
		Chamber_x[i] = i;
		Chamber_y[i] = Dc/2;
	}
	for(i = (int)x[1]+1; i<x[2];i++){
		Chamber_x[i] = i;
		Chamber_y[i] = pow(pow(size1,2)-pow(i-x[1],2),0.5) + Dc/2 - size1;
	}
	for(i = (int)x[2]+1; i<x[3];i++){
		Chamber_x[i] = i;
		Chamber_y[i] = (y[3]-y[2])/(x[3]-x[2])*i + y[2] - (y[3]-y[2])/(x[3]-x[2])*x[2];
	}
	for(i = (int)x[3]+1; i<x[4];i++){
		Chamber_x[i] = i;
		Chamber_y[i] = -pow(pow(Dt / 2 * 1.5,2)-pow(i-x[4],2),0.5) + Dt/2 + Dt / 2 * 1.5;
	}
	

//	for(i = 0; i < x[6];i++){
//		printf("x = %f\ty = %f\n",Chamber_x[i],Chamber_y[i]);
//	}
	
	
	fp = fopen("chamber.dxf","w");
	dxf_header(&fp);
	dxf_line(&fp,x[0],y[0],x[1],y[1]);
	dxf_arc(&fp,x[1],y[1]-size1,size1,90-angle1,90);
	dxf_line(&fp,x[2],y[2],x[3],y[3]);
	dxf_arc(&fp,x[4],y[4]+Dt/2*1.5,Dt/2*1.5,270-angle1,270);
	fflush(fp);
	fclose(fp);

	return x[6];
}

void calc_conical_nozzle(void){

	int i;
	FILE *fp;
	double x[3],y[3];
	int exit_angle = 15;

	x[0] = throat_axis;
	y[0] = Dt / 2;
	x[1] = Dt / 2 * 1.5 * cos(convert_to(rad,270 + exit_angle)) + x[0];
	y[1] = Dt / 2 * 1.5 * sin(convert_to(rad,270 + exit_angle)) + Dt / 2 + Dt / 2 * 1.5;
	x[2] = x[0] + Ln;
	y[2] = (x[2]-x[1])*tan(convert_to(rad,exit_angle)) + y[1];

	for(i = (int)x[0]; i<x[1];i++){
		Chamber_x[i] = i;
		Chamber_y[i] = -pow(pow(Dt / 2 * 1.5,2)-pow(i-x[0],2),0.5) + Dt/2 + Dt / 2 * 1.5;
	}
	for(i = (int)x[1]+1;i <x[2];i++){
		Chamber_x[i] = i;
		Chamber_y[i] = (y[2]-y[1])/(x[2]-x[1])*i + y[1] - (y[2]-y[1])/(x[2]-x[1])*x[1];
	}
	exit_axis = (int)x[2];

	fp = fopen("chamber.dxf","a");
	dxf_arc(&fp,x[0],y[0]+Dt/2*1.5,Dt/2*1.5,270,270+exit_angle);
	dxf_line(&fp,x[1],y[1],x[2],y[2]);
	dxf_footer(&fp);
	fflush(fp);
	fclose(fp);
}

void calc_foelsch_nozzle(void){
	int i,j;
	double theta;
	double tmp_theta;
	double ve,v1,mach,A,gamma;
	double y0,y1,x1,r,r0,r1,alpha;
	double l,l1,l2,l3;
	double angle;

	theta = convert_to(rad,12);
	
	y0 = Dt/2;
	ve = prandtle_meyer(CEA[nozzle_exit].Mach);
	v1 = 2 * (ve / 2 - theta);
	gamma = CEA[chamber].Gamma;

	mach = get_mach_from_prandtle_meyer(v1);

	A = (1/mach)*pow(((gamma-1)*mach*mach + 2)/(gamma + 1),((gamma + 1)/(2*(gamma -1))));
	r1 = y0*sqrt(A/(2*(1-cos(theta))));
	y1 = r1*sin(theta);
	x1 = 3*0.5*(y1-y0)/tan(theta);

	j = throat_axis;
	for(i = 0; i<x1; i++){
		Chamber_x[j] = i + throat_axis;
		Chamber_y[j] = y0 + (tan(theta)/x1)*pow(i,2)*(1-i/(3*x1));
//		printf("x = %f, y =%f\n",Chamber_x[j],Chamber_y[j]);
		j++;
	}

	r0 = r1/sqrt(A);
//	printf("A = %f\n",A);
//	printf("y1 = %f\n",y1);
//	printf("x1 = %f\n",x1);
//	printf("r0 = %f\n",r0);
//	printf("r1 = %f\n",r1);

// from here, increment mach number 
	for(;CEA[nozzle_exit].Mach > mach; mach = mach + 0.005){
		r = r0*sqrt((1/mach)*
			pow(((gamma-1)*mach*mach + 2)/(gamma + 1),((gamma + 1)/(2*(gamma -1)))));
		alpha = asin(1/mach);
		tmp_theta = 0.5*(ve-prandtle_meyer(mach));
		l1 = 2*r*sin(alpha)*sin(tmp_theta);
		l2 = sqrt(pow(2*r*sin(alpha)*sin(tmp_theta),2) + 
			8*r*r*sin(alpha)*(cos(tmp_theta)-cos(theta)) * 
			(cos(alpha)*sin(tmp_theta) + sin(alpha + tmp_theta)));
		l3 = (-2*sin(alpha)*(cos(alpha)*sin(tmp_theta) + sin(alpha + tmp_theta)));
		l = (l1 - l2)/(l3);
	
		Chamber_x[j] = r*cos(tmp_theta) + l*cos(tmp_theta + alpha) +
				(x1-r1*cos(theta)) + throat_axis;
		Chamber_y[j] = l*sin(tmp_theta + alpha) + r*sin(tmp_theta);
		
		
//		printf("r = %f\n",r);
//		printf("alpha = %f\n",alpha);
//		printf("tmp_theta = %f\n",tmp_theta);
//		printf("l = %f\n",l);
//		printf("l1 = %f\n",l1);
//		printf("l2 = %f\n",l2);
//		printf("l3 = %f\n",l3);
//		printf("x = %f, y = %f, j = %d\n",Chamber_x[j],Chamber_y[j],j);
		j++;

	}
	angle = atan2((Chamber_y[j-1]-Chamber_y[j-2]),(Chamber_x[j-1]-Chamber_x[j-2]));	
	exit_axis = j-1;
	printf("angle = %f\n",convert_to(deg,angle));
	printf("exit = %d\n",exit_axis);

}

//can calculate from mach 1-10
double get_mach_from_prandtle_meyer(double v1){

	double num;
	double mach;

	mach = 1;
	for(num = 0; num < 7; num++){	//change here for mach precision
		while(1){
			printf("v1 = %f < pr = %f mach = %f\n",v1,prandtle_meyer(mach),mach);
			if(v1<prandtle_meyer(mach))break;
			mach = mach + 1/pow(10,num);
		}
		mach = mach - 1/pow(10,num);
	}
	printf("mach = %f\n",mach);
	return mach;
}

double prandtle_meyer(double mach){
	double vm;
	double gamma;

	gamma = CEA[chamber].Gamma;

	vm = sqrt((gamma+1)/(gamma-1)) * 
		atan(sqrt(((gamma-1)/(gamma+1))*(mach*mach - 1))) - 
		atan(sqrt(mach*mach - 1));

	return vm;
}

void calc_rao_nozzle_(void){
}



void calc_nozzle_geom(void){
	
}


double CEA_coeff(int x, double chamber, double throat, double nozzle_exit){

	double diff_c1;
	double coeff;
	double val;

	if(x < throat_axis){
		diff_c1 = (Dc - Dt)/2;
		coeff = (Chamber_y[x] - Dt/2)/diff_c1;
		val = (chamber - throat)*coeff + throat;
	}else if(x >= throat_axis){
		diff_c1 = (De - Dt)/2;
		coeff = (Chamber_y[x] - Dt/2)/diff_c1;
		val = (nozzle_exit - throat)*coeff + throat;
	}

//	printf("coeff = %f\n",val);

	return val;
}

void plot_chamber(void){
	
	FILE *gp;
	int i;

	gp = popen("gnuplot -persist","w");
	fprintf(gp,"set datafile separator \"\t\"\n");
	fprintf(gp,"set size square\n");
	fprintf(gp,"set title 'Kerosene'\n");
	fprintf(gp,"set xlabel '[mm]'\n");
	fprintf(gp,"set ylabel '[mm]'\n");
	fprintf(gp,"set xrange[-10:%f]\n",Chamber_x[exit_axis]+10);
	fprintf(gp,"set yrange[-10:%f]\n",Chamber_x[exit_axis]+10);
	fprintf(gp,"plot '-' with lines\n");
	
	for(i = 0; i<exit_axis; i++){
		fprintf(gp,"%f\t%f\n",Chamber_x[i],Chamber_y[i]);
	}
	fprintf(gp,"e\n");

	pclose(gp);
}

