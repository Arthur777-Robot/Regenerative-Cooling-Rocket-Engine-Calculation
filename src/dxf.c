#include<stdio.h>
#include"dxf.h"

void dxf_header(FILE **fp){

	fprintf(*fp,"0\nSECTION\n"
			"2\nTABLES\n"
			"0\nTABLE\n"
			"0\nLAYER\n"
			"2\nLayer0\n"
			"70\n1\n"
			"62\n7\n"
			"6\nCONTINUOUS\n"
			"0\nENDTAB\n"
			"0\nENDSEC\n"
			"0\nSECTION\n"
			"2\nENTITIES\n");

}

void dxf_footer(FILE **fp){
	fprintf(*fp, "0\n"
			"ENDSEC\n"
			"0\n"
			"EOF\n");
}

void dxf_line(FILE **fp, float x1, float y1, float x2, float y2){
	fprintf(*fp,"0\n"
			"LINE\n"
			"8\n"
			"Layer0\n"
			"6\n"
			"CONTINUOUS\n"
			"62\n"
			"8\n");
	fprintf(*fp, "10\n%f\n20\n%f\n11\n%f\n21\n%f\n",x1,y1,x2,y2);
}

void dxf_arc(FILE **fp, float x, float y, float r, float theta1, float theta2){
	fprintf(*fp,"0\n"
			"ARC\n"
			"8\n"
			"Layer0\n"
			"6\n"
			"CONTINUOUS\n"
			"62\n"
			"8\n");
	fprintf(*fp, "10\n%f\n20\n%f\n40\n%f\n50\n%f\n51\n%f\n",x,y,r,theta1,theta2);
}

void dxf_circle(FILE **fp, float x, float y, float r){
	fprintf(*fp,"0\n"
			"CIRCLE\n"
			"8\n"
			"Layer0\n"
			"6\n"
			"CONTINUOUS\n"
			"62\n"
			"8\n");
	fprintf(*fp, "10\n%f\n20\n%f\n40\n%f\n",x,y,r);
}
