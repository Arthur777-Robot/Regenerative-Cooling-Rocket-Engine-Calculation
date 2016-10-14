#include <stdio.h>
#include <math.h>
#include "valuables.h"

float convert_to(int prefix, float num){

	if((prefix == mm) || (prefix == gram)){
		return num * 1000;

	}else if((prefix == m) || (prefix == kg) ||(prefix == L)){
		return num / 1000;

	}else if((prefix == mm2) ||(prefix == Pa)){
		return num * 1000000;

	}else if((prefix == m2) || (prefix == MPa)){
		return num / 1000000;

	}else if(prefix == mm3){
		return num * 1000000000;

	}else if(prefix == m3){
		return num / 1000000000;

	}else if(prefix == deg_c){
		return num - 273.15;
	
	}else if(prefix == deg_k){
		return num + 273.15;
	
	}else if(prefix == rad){
		return num * M_PI / 180;
	
	}else if(prefix == deg){
		return num / M_PI * 180;
	
	}else{
		return 0;
	}
}


