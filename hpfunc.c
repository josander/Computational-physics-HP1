/*
hpfunc.c
Contains functions for homeproblem 1/b

By Jossan och Svensson

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "alpotential.h"
#define PI 3.141592653589
#define K_B 0.000086173324
#define nbr_of_atoms 256

// Randomly displaces the atoms in position[][3] by 5% of the lattice parameter
void rand_disp(double position[][3] ,double lattice_param ,int N)
{
	int i, j;
	double sum;	
	double ran[3];
	double length_corr;	
	for(i = 0; i < N; i++){
		sum = 0.0;		
		for(j = 0; j < 3; j++){
			ran[j] = ((double) rand() / (double) RAND_MAX) -0.5; // ran[j] goes from -0.5 to 0.5
			sum += ran[j]*ran[j];
		}
		
		sum = sqrt(sum);
	
		length_corr = (lattice_param*0.05)/sum;
				
		
		for(j = 0; j < 3; j++){
			position[i][j] +=  ran[j] * length_corr;
		}
	}  		  
}

// Calculates kinetic energy of the system
double get_ke(double v[][3], int N, double m)
{
	int i, j;
	double ke = 0.0;

	for (i = 0; i < N; i++){
		for(j = 0; j < 3; j++){		
			ke += 0.5*m*v[j][i]*v[j][i];
		}
	}
	return ke;
}

// Calculates temperature
double get_T(double ke, int N)
{
	
	double T = 2*ke/((3*N-5)*K_B);
	return T;	
}

// Calculates alpha_T, the correction parameter for v_i
void rescale_T(double timestep, double tau_T, double T_eq, double T, double v[][3], int N)
{
	double alphaT;
	int i, j;

	alphaT = 1 + timestep*(T_eq - T)/(tau_T*T);

	for(i = 0; i < N; i++){
		for(j = 0; j < 3; j++){
			v[i][j] = sqrt(alphaT) * v[i][j];
		}
	}
}

// Calculates the preassure
double get_P(double q[][3], double cell_length, int N, double T)
{
	double virial, P;	
	virial = get_virial_AL(q, cell_length, N);
	
	P = (N*K_B*T + virial)/pow(cell_length,3);

	return P;

}
// Rescales lattice_param according to eqlibr
double rescale_P(double timestep, double tau_P, double P_eq, double P, double q[][3], int N, double kappa_T, double lattice_param)
{
	double alphaP;
	int i, j;


	alphaP = pow((1.0 - timestep*kappa_T*(P_eq - P)/tau_P), 0.333333);	

	lattice_param =  alphaP*lattice_param;
	//printf("%F \n", alphaP);
	for(i = 0; i < N; i++){
		for(j = 0; j < 3; j++){
			q[i][j] = alphaP * q[i][j];
		}
	}
	return lattice_param;
}


// Function that calculates the correlation function, the statistical inefficiency and the standard deviation
void get_corr_func(double A[], double *corr_func, int nbr_of_timesteps, int start)
{
	int i, k;
	double mean = 0;
	double mean2 = 0;
	int stop = 1000;
	double first_term[nbr_of_timesteps-start];
	double s = 0;
	double sigmaTot;

	// Initiate the array first_term
	for(k = 0; k < (nbr_of_timesteps - start); k++){
		first_term[k] = 0.0;
	}

	// Calculate all the expected values of A
	for(i = start; i < nbr_of_timesteps; i++){
		mean += A[i]/(nbr_of_timesteps - start);
		mean2 += ((A[i]*A[i])/(nbr_of_timesteps - start)); 
	}

	// Calculate the first term
	for(i = start; i < nbr_of_timesteps; i++){
		for(k = 0; k < (nbr_of_timesteps-i); k++){
			first_term[k] += (A[i]*A[i+k])/(nbr_of_timesteps - start - k);
		}
	}

	// Calculate the correlation function
	for(k = 0; k < (nbr_of_timesteps - start); k++){
		corr_func[k] = ((first_term[k] - (mean*mean))/(mean2 - (mean*mean)));
	}

	// Calculate the statistical inefficiency
	i = 0;
	while(corr_func[i] >= exp(-2)){
		s += corr_func[i];
		i++;
	}

	s *= 2;

	sigmaTot = sqrt((mean2 - mean*mean)/(nbr_of_timesteps-start)*s);

	printf("Statistical inefficiency: %F \n", s);

	printf("Result: %.8e ± %.10e \n", mean, sigmaTot);


}

// Function that calculates the spectral function
void get_spectral_func(double *vel_corr_func, double *omega, double *spectral_func, int nbr_of_steps, int num_of_freq, int timestep)
{
	int i, j;
	double factor = 2.0*5/nbr_of_steps; // 5 = nbr_of_step*timestep
	for(i = 0; i < num_of_freq ; i++){
		omega[i] = i * PI / num_of_freq;
		spectral_func[i] = 0;		
		for(j = 0; j < nbr_of_steps; j++){
			spectral_func[i] += vel_corr_func[j] * cos(omega[i] * j);
		}
		spectral_func[i] *= factor;
	}
}

 












