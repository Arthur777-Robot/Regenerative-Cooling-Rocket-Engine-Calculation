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

//equilibrium and frozen option must be written
void get_CEA_param(void){

	char temp1[256],temp2[1];
	FILE *fp;

	char s[8][10] = {"T, K","Cp, KJ","VISC","PRANDTL","Ae/At","CSTAR","CF","Isp"};
	int i = 0;
	char split[] = " \n";
	char *addr1,*addr2;

	while(1){
		printf("calculate with frozen or equilibrium option\n if frozen, enter f. if equilibrium, enter e\n");
		scanf("%c",temp2);

		if(temp2[0] == 'f'){
			printf("calculate with frozen option\n");
			break;
		}else if(temp2[0] == 'e'){
			printf("calculate with equilibrium option\n");
			break;
		}

	}

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
					printf("Tc\n");
					CEA[0].Tc = atof(addr2);
					printf("%f\n",CEA[0].Tc);
					CEA[1].Tc = atof(strtok(NULL,split));
					printf("%f\n",CEA[1].Tc);
					CEA[2].Tc = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].Tc);
					break;
				}else if(i == 1){
					printf("Cp\n");
					CEA[0].Cp = atof(addr2);
					printf("%f\n",CEA[0].Cp);
					CEA[1].Cp = atof(strtok(NULL,split));
					printf("%f\n",CEA[1].Cp);
					CEA[2].Cp = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].Cp);
					i--;
					break;
				}else if(i == 2){
					printf("Visc_gas\n");
					CEA[0].Visc_gas = atof(addr2);
					printf("%f\n",CEA[0].Visc_gas);
					CEA[1].Tc = atof(strtok(NULL,split));
					printf("%f\n",CEA[1].Visc_gas);
					CEA[2].Tc = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].Visc_gas);
					break;
				}else if(i == 3){
					printf("Prandtl_gas\n");
					CEA[0].Prandtl_gas = atof(addr2);
					printf("%f\n",CEA[0].Prandtl_gas);
					CEA[1].Prandtl_gas = atof(strtok(NULL,split));
					printf("%f\n",CEA[1].Prandtl_gas);
					CEA[2].Prandtl_gas = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].Prandtl_gas);
					break;
				}else if(i == 4){
					printf("AeAt\n");
					CEA[0].AeAt = 0;
					printf("%f\n",CEA[0].AeAt);
					CEA[1].AeAt = atof(addr2);
					printf("%f\n",CEA[1].AeAt);
					CEA[2].AeAt = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].AeAt);
					break;
				}else if(i == 5){
					printf("Cstar\n");
					CEA[0].Cstar = 0;
					printf("%f\n",CEA[0].Cstar);
					CEA[1].Cstar = atof(addr2);
					printf("%f\n",CEA[1].Cstar);
					CEA[2].Cstar = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].Cstar);
					break;
				}else if(i == 6){
					printf("Cf\n");
					CEA[0].Cf = 0;
					printf("%f\n",CEA[0].Cf);
					CEA[1].Cf = atof(addr2);
					printf("%f\n",CEA[1].Cf);
					CEA[2].Cf = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].Cf);
					break;
				}else if(i == 7){
					printf("ISP\n");
					CEA[0].Isp = 0;
					printf("%f\n",CEA[0].Isp);
					CEA[1].Isp = atof(addr2);
					printf("%f\n",CEA[1].Isp);
					CEA[2].Isp = atof(strtok(NULL,split));
					printf("%f\n",CEA[2].Isp);
					break;
				}
			}
			if(i > 7) break;
			i++;
		}
	}
}

