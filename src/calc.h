#ifndef CALC_H
#define CALC_H

void calc_chamber_spec(void);
void calc_nozzle(void);
void calc_chamber(void);
int calc_chamber_geom(void);
void calc_nozzle_geom(void);
void calc_chamber_strength(void);
void calc_fuel_consumption(void);
void calc_rao_nozzle(void);
void calc_conical(void);
void calc_foelsch_nozzle(void);
float prandtle_meyer(float);
float specific_heat_ratio(void);

// using bartz
void calc_gas_heat_transfer(float,float,float,float);
void calc_fuel_heat_transfer(void);
void calc_metal_heat_transfer(void);
void calc_total_heat_transfer(int,float, float*, float*, float*);
float calc_delta_fuel_temp(float, float);
void calc_fuel_cost(void);

float CEA_coeff(int, float, float, float);
void plot_chamber(void);

#endif
