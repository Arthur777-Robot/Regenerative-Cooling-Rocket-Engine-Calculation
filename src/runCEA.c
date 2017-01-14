#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "runCEA.h"
#include "valuables.h"

T_CEA CEA[3]; // parameter for chamber, throat, exit

void set_CEA(void){
	make_inp();
	run_CEA();
	get_CEA_param();
}

// rewrite CEA input file here to your required parameters
void make_inp(void){

	FILE *fp;

	fp = fopen("fuel.inp","w");
	fprintf(fp,"prob\n"
				"rocket equilibrium frozen nfz=2\n"
				"o/f %f\n"
				"p,bar %f\n"
				"pip %f\n"
				"reac\n"
				" fuel %s wt%%=100 t,k=298.15\n"
				" oxid %s wt%%=100 t,k=90.17\n"
				"output transport\n"
				"output short\n"
				"output trace=1e-5\n"
				"plot p pi/p o/f isp ivac aeat cf t\n"
				"end\n",
				(float)OF,(float)(PC*10),(float)(PC*10),
				TYPE_FUEL,TYPE_OXIDIZER);

	fflush(fp);
	fclose(fp);
	
}

void run_CEA(void){
	FILE *fp;
	fp = popen("FCEA2 > /dev/null","w");
	fprintf(fp,"fuel");
	pclose(fp);
}

//hate this code. please someone rewite this with smart brain.
void get_CEA_param(void){

	char temp1[256],temp2[1];
	FILE *fp;

	char s[11][12] = {"T, K","MACH NUMBER","VISC","Cp, KJ","PRANDTL","Cp, KJ","PRANDTL","Ae/At","CSTAR","CF","Isp"};
	int i = 0, j = 0;
	char split[] = " \n";
	char *addr1,*addr2;

	float cp1[3],cp2[3];
	float pr1[3],pr2[3];

	if((fp = fopen("fuel.out","r")) == NULL){
		printf("error reading fuel.out file\n");
	}

	while(fgets(temp1,256,fp)!=NULL){
		if((strstr(temp1,s[i]))!= NULL){
			addr1 = strstr(temp1,"  ");
			addr2 = strtok(addr1,split);
			while(addr2 != NULL){
//				printf("%s\n",addr2);
//				addr2 = strtok(NULL,split);
				if(i == 0){
					CEA[0].Tc = atof(addr2);
					CEA[1].Tc = atof(strtok(NULL,split));
					CEA[2].Tc = atof(strtok(NULL,split));

//					printf("Tc\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",CEA[j].Tc);
//					}
					break;
				}else if(i == 1){
					CEA[0].Mach = atof(addr2);
					CEA[1].Mach = atof(strtok(NULL,split));
					CEA[2].Mach = atof(strtok(NULL,split));
//					printf("Mach\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",CEA[j].Mach);
//					}
					break;

				}else if(i == 2){
					CEA[0].Visc_gas = atof(addr2);
					CEA[1].Visc_gas = atof(strtok(NULL,split));
					CEA[2].Visc_gas = atof(strtok(NULL,split));

//					printf("Visc_gas\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",CEA[j].Visc_gas);
//					}
					break;
				}else if(i == 3){
					cp1[0]= atof(addr2);
					cp1[1]= atof(strtok(NULL,split));
					cp1[2] = atof(strtok(NULL,split));

//					printf("Equilibrium Cp\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",cp1[j]);
//					}
					break;
				}else if(i == 4){
					pr1[0] = atof(addr2);
					pr1[1] = atof(strtok(NULL,split));
					pr1[2] = atof(strtok(NULL,split));

//					printf("Equilibrium Prandtl_gas\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",pr1[j]);
//					}
					break;
				}else if(i == 5){
					cp2[0] = atof(addr2);
					cp2[1] = atof(strtok(NULL,split));
					cp2[2] = atof(strtok(NULL,split));

//					printf("Frozen Cp\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",cp2[j]);
//					}
					break;
				}else if(i == 6){
					pr2[0] = atof(addr2);
					pr2[1] = atof(strtok(NULL,split));
					pr2[2] = atof(strtok(NULL,split));

//					printf("Frozen Prandtl_gas\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",pr2[j]);
//					}
					break;
				}else if(i == 7){
					CEA[0].AeAt = 0;
					CEA[1].AeAt = atof(addr2);
					CEA[2].AeAt = atof(strtok(NULL,split));

//					printf("AeAt\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",CEA[j].AeAt);
//					}
					break;
				}else if(i == 8){
					CEA[0].Cstar = 0;
					CEA[1].Cstar = atof(addr2);
					CEA[2].Cstar = atof(strtok(NULL,split));

//					printf("Cstar\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",CEA[j].Cstar);
//					}
					break;
				}else if(i == 9){
					CEA[0].Cf = 0;
					CEA[1].Cf = atof(addr2);
					CEA[2].Cf = atof(strtok(NULL,split));

//					printf("Cf\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",CEA[j].Cf);
//					}
					break;
				}else if(i == 10){
					CEA[0].Isp = 0;
					CEA[1].Isp = atof(addr2);
					CEA[2].Isp = atof(strtok(NULL,split));

//					printf("ISP\n");
//					for(j = 0;j < 3; j++){
//						printf("%f\n",CEA[j].Isp);
//					}
					break;
				}
			}
			if(i > 10) break;
			i++;
		}
	}
	while(1){
		printf("calculate with frozen or equilibrium option\n if frozen, enter f. if equilibrium, enter e\n");
		scanf("%c",temp2);

		if(temp2[0] == 'f'){
			printf("calculate with frozen option\n");
			for(i = 0; i < 3; i++){
				CEA[i].Cp = cp2[i];
				CEA[i].Prandtl_gas = pr2[i];
			}
			break;
		}else if(temp2[0] == 'e'){
			printf("calculate with equilibrium option\n");
			for(i = 0; i < 3; i++){
				CEA[i].Cp = cp1[i];
				CEA[i].Prandtl_gas = pr1[i];
			}
			break;
		}else{
			printf("error, please type f or e");
		}

	}
}

