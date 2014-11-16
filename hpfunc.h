/*
hpfunc.h
 
.
*/

#ifndef _hpfunc_h
#define _hpfunc_h

//extern void get_forces_AL(double[][3] , double[][3], double, int);
//extern double get_energy_AL(double[][3], double, int);
//extern double get_virial_AL(double[][3], double, int);
extern void rand_disp(double[][3],double, int);
extern double get_ke(double[][3], int, double);
extern double get_T(double, int);
extern double get_P(double[][3], double ,int , double)
extern void rescale_T(double, double, double , double, double[][3], int);
void rescale_P(double, double, double, double, double , int , double )


#endif
