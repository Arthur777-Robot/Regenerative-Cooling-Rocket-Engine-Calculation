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
void calc_conical_nozzle(void);
void calc_foelsch_nozzle(void);
double get_mach_from_prandtle_meyer(double);
double prandtle_meyer(double);

// using bartz
void calc_gas_heat_transfer(int,double,double,double);
void calc_fuel_heat_transfer(void);
void calc_metal_heat_transfer(void);
void calc_total_heat_transfer(int,double, double*, double*, double*);
double calc_delta_fuel_temp(double, double);
void calc_fuel_cost(void);

double CEA_coeff(int, double, double, double);
void boundary_layer(int);
void plot_chamber(void);
double rp1_rho(double, double);
double rp1_visc(double);
double rp1_cp(double);

#endif
