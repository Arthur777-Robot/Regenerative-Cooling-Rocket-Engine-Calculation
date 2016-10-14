#ifndef CALC_H
#define CALC_H

void calc_chamber_spec(void);
void calc_nozzle(void);
void calc_chamber(void);
void calc_chamber_strength(void);
void calc_fuel_consumption(void);
void calc_regene(void);

// using bartz
void calc_gas_heat_stansfer_coeff(void);
void calc_fuel_heat_stansfer_coeff(void);
void calc_metal_heat_stansfer_coeff(void);
void calc_total_heat_stansfer_coeff(void);
void calc_delta_fuel_temp(void);
void calc_fuel_cost(void);

#endif
