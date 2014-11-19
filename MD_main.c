/*
 MD_main.c
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"
#include "hpfunc.h"
#define PI 3.141592653589
#define K_B 0.000086173324

/* OBS: Kolla om vi får tillbaka väntevärdet av virialen eller ej! */


// Main program 
int main()
{

    	// Declaration of variables 
	double lattice_param, cell_size;
	double timestep;
	int nbr_of_timesteps;
	int nbr_of_atoms;
	int Nx, Ny, Nz;
     	double random_value;
	int i, j, n;
	double m;
	double energy, pe, ke;
	double tau_T, tau_P;
	double temp_eq;
	double press_eq;
	double kappa_P;
	double s_T;
	double s_P;
	int startCut;


	// Initiation of variables 
	lattice_param = 4.05; // Units: [Å]
	timestep = 0.01; // [ps]
	nbr_of_timesteps = 2000;
	nbr_of_atoms = 256;
	Nx = 4, Ny = 4, Nz = 4;
	m = 0.00279636665; // Metal units [ev/Å]
	temp_eq = 500 + 273.15; // Degree Celsius 
	press_eq = 6.324209 * pow(10, -7); // 1 Atm in eV/Å^3
	tau_T = timestep*100;
	tau_P = timestep*100;
	kappa_P = 2.21901454; //3.85 * pow(10, 9);/ // Liquid Aluminum Units: Å^3/eV
	cell_size = lattice_param*Nx;
	s_T = 0; 
	s_P = 0;
	startCut = 800;

	// Declaration of matrixes and arrays 
	double q[4*Nx*Ny*Nz][3], v[nbr_of_atoms][3], a[nbr_of_atoms][3];
	double f[4*Nx*Ny*Nz][3];
	double *temp = malloc(nbr_of_timesteps * sizeof(double));
	double *press = malloc(nbr_of_timesteps * sizeof(double));
	double *corr_func_T = malloc((nbr_of_timesteps-startCut+1) * sizeof(double));
	double *corr_func_P = malloc((nbr_of_timesteps-startCut+1) * sizeof(double));

	// Initiation of corr_func
	for(i = 0; i <nbr_of_timesteps; i++){
		corr_func_T[i] = 0.0;
		corr_func_P[i] = 0.0;
	}

	// Initiation of the fcc lattice of Al-atoms 
	init_fcc(q, Nx, lattice_param);

	// Introduction of displacements, about 5 % of the lattice spacing
     	srand(time(NULL)); // Should only be called once
	rand_disp(q, lattice_param, nbr_of_atoms);

	// Initiation of the velocities 
	for(j = 0; j < nbr_of_atoms; j++){
		for(n = 0; n < 3; n++){
			v[j][n] = 0.0;
		}
	}

	// Get forces
	get_forces_AL(f, q, cell_size, nbr_of_atoms);

	// Scale forces to acceleration
	for(j = 0; j < nbr_of_atoms; j++){
		for(n = 0; n < 3; n++){
			a[j][n] = f[j][n]/m;
		}
	}

	// Calculate of initial energies
	pe = get_energy_AL(q, cell_size, nbr_of_atoms);
	pe = sqrt(pe*pe);
	ke = get_ke(v, nbr_of_atoms, m);
	energy = sqrt(pe*pe) + sqrt(ke*ke);

	// Calculate initial temperature
	temp[0] = get_T(ke, nbr_of_atoms);

	// Calculate initial pressure
	press[0] = get_P(q, cell_size, nbr_of_atoms, temp[0]);

	// Make a file to save the energies in
	FILE *e_file;
	e_file = fopen("energy.data","w");

	// Save the initial energies in the file
	fprintf(e_file,"%.5f \t %e \t %e \t %e \t %F \t %e \n", 0.0, energy, pe, ke, temp[0], press[0]);

	// Time evolution according to the velocity Verlet algorithm
	for (i = 1; i < nbr_of_timesteps + 1; i++){
		// v(t+dt/2)
		for (j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
		    		v[j][n] += timestep * 0.5 * a[j][n];
			}
		} 

		// q(t+dt) 
		for (j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
		    		q[j][n] += timestep * v[j][n];
			}
		}

		// a(t+dt)
		// Get forces
		get_forces_AL(f, q, cell_size, nbr_of_atoms);

		// Scale forces to acceleration
		for(j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
				a[j][n] = f[j][n]/m;
			}
		}

		// v(t+dt)
		for (j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
		    		v[j][n] += timestep * 0.5 * a[j][n];
			}
		} 

		// Calculate ke to calculate rescale
		ke = get_ke(v, nbr_of_atoms, m);

		// Calculate the temperature
		temp[i] = get_T(ke, nbr_of_atoms);

		// Scale velocity of the atoms to obtain the right temperature
		rescale_T(timestep, tau_T, temp_eq, temp[i], v, nbr_of_atoms);

		// Calculate the pressure
		cell_size = lattice_param * Nx;
		press[i] = get_P(q, cell_size, nbr_of_atoms, temp[i]);
	
		// Scale position of the atoms to obtain the right pressure
		lattice_param = rescale_P(timestep, tau_P, press_eq, press[i], q, nbr_of_atoms, kappa_P, lattice_param);

		// Calcutaion of the pe, ke and total energy
		pe = get_energy_AL(q, cell_size, nbr_of_atoms);
		pe = sqrt(pe*pe);
		ke = get_ke(v, nbr_of_atoms, m);
		energy = sqrt(pe*pe) + sqrt(ke*ke);

		// Print every 1000 timestep in the terminal
		if(i%1000 == 0){
			printf("%i av %i steg \n", i, nbr_of_timesteps);
		}
	
		// Print the average energy data to output file
		fprintf(e_file,"%.5f \t %e \t %e \t %e \t %F \t %e \n", i*timestep, energy, pe, ke, temp[i], press[i]);
	}

	// Get the correlation functions
	get_corr_func(temp, corr_func_T, nbr_of_timesteps+1, startCut);
	get_corr_func(press, corr_func_P, nbr_of_timesteps+1, startCut);

	// Write corr_func to data file
	FILE *c_file;
	c_file = fopen("correlation.data","w");

	for(i = 0; i < (nbr_of_timesteps-startCut+1); i++){
		fprintf(c_file,"%.5f \t %e \n", corr_func_T[i], corr_func_P[i]);
	}

	// Calculate the statistical inefficiency
	for(i = 0; i < 500; i++){
		s_T += corr_func_T[i];
	}
	for(i = 0; i < 700; i++){
		s_P += corr_func_P[i];
	}

	s_T *= 2;
	s_P *= 2;

	printf("sT: %F \t sP: %F \n", s_T, s_P);

	// Close the energy output file 
	fclose(e_file);
	fclose(c_file);

	// Free allocated memory
	free(temp); free(press); free(corr_func_T); free(corr_func_P); 
	temp = NULL; press = NULL; corr_func_T = NULL; corr_func_P = NULL;

	return 0;
}
