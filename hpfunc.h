/*
hpfunc.h
 
.
*/

#ifndef _hpfunc_h
#define _hpfunc_h

extern void rand_disp(double[][3],double, int);
extern double get_ke(double[][3], int, double);
extern double get_T(double, int);
extern void rescale_T(double, double, double , double, double[][3], int);
extern double get_P(double[][3], double ,int , double);
extern double rescale_P(double, double, double, double, double[][3], int , double, double );
extern void get_corr_func(double[], double *, int);


#endif
