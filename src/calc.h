#ifndef CALC_H
#define CALC_H

void calc_chamber_spec(void);
void calc_nozzle(void);
void calc_chamber(int*,float*);
void calc_chamber_strength(void);
void calc_fuel_consumption(void);

// using bartz
void calc_gas_heat_transfer(float,float,float,float);
void calc_fuel_heat_transfer(void);
void calc_metal_heat_transfer(void);
void calc_total_heat_transfer(int,float, float*, float*, float*);
float calc_delta_fuel_temp(float, float);
void calc_fuel_cost(void);

int plot_chamber(int*, float*);
float CEA_coeff(int,float*, float, float, float);

#endif
