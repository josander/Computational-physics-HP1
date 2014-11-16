/*
 MD_main.c
 
 Created by AL on 2013-10-31.
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


// Main program 
int main()
{

    	// Declaration of variables 
	double lattice_param;
	double timestep;
	int nbr_of_timesteps;
	int nbr_of_atoms;
	int Nx, Ny, Nz;
     	double random_value;
	int i, j, n;
	double m;
	double energy, pe, ke;
	double tau_T, tau_P;
	double temp;
	double temp_eq;
	double press;
	double press_eq;
	double kappa_P;


	// Initiation of variables 
	lattice_param = 4.05; // Units: [Å]
	timestep = 0.01; // [ps]
	nbr_of_timesteps = 15000;
	nbr_of_atoms = 256;
	Nx = 4, Ny = 4, Nz = 4;
	m = 0.00279636665; // Metal units [ev/Å]
	temp_eq = 500 + 273.15; // Degree Celsius
	press_eq = 6.324209 * pow(10, -7); // 1 Atm in eV/Å^3
	tau_T = timestep*100;
	tau_P = timestep*2000;
	kappa_P = 2.21901454; //3.85 * pow(10, 9); // Liquid Aluminum Units: Å^3/eV

	// Declaration of matrixes and arrays 
	double q[4*Nx*Ny*Nz][3], v[nbr_of_atoms][3], a[nbr_of_atoms][3];
	double f[4*Nx*Ny*Nz][3];

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
	get_forces_AL(f, q, Nx*lattice_param, nbr_of_atoms);

	// Scale forces to acceleration
	for(j = 0; j < nbr_of_atoms; j++){
		for(n = 0; n < 3; n++){
			a[j][n] = f[j][n]/m;
		}
	}

	// Calculate of initial energies
	pe = get_energy_AL(q, Nx*lattice_param, nbr_of_atoms);
	pe = sqrt(pe*pe);
	ke = get_ke(v, nbr_of_atoms, m);
	energy = sqrt(pe*pe) + sqrt(ke*ke);

	// Calculate initial temperature
	temp = get_T(ke, nbr_of_atoms);

	// Calculate initial pressure
	press = get_P(q, Nx*lattice_param, nbr_of_atoms, temp);

	// Make a file to save the energies in
	FILE *e_file;
	e_file = fopen("energy.data","w");

	// Save the initial energies in the file
	fprintf(e_file,"%.5f \t %e \t %e \t %e \t %F \t %e \n", 0.0, energy, pe, ke, temp, press);

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
		get_forces_AL(f, q, Nx*lattice_param, nbr_of_atoms);

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
		temp = get_T(ke, nbr_of_atoms);

		// Scale velocity of the atoms to obtain the right temperature
		rescale_T(timestep, tau_T, temp_eq, temp, v, nbr_of_atoms);

		// Calculate the pressure
		press = get_P(q, Nx*lattice_param, nbr_of_atoms, temp);

		// Scale position of the atoms to obtain the right pressure
		rescale_P(timestep, tau_P, press_eq, press, q, nbr_of_atoms, kappa_P);

		// Calcutaion of the pe, ke and total energy
		pe = get_energy_AL(q, Nx*lattice_param, nbr_of_atoms);
		pe = sqrt(pe*pe);
		ke = get_ke(v, nbr_of_atoms, m);
		energy = sqrt(pe*pe) + sqrt(ke*ke);

		// Print every 1000 timestep in the terminal
		if(i%1000 == 0){
			printf("%i av %i steg \n", i, nbr_of_timesteps);
		}
	
		// Print the average energy data to output file
		fprintf(e_file,"%.5f \t %e \t %e \t %e \t %F \t %e \n", i*timestep, energy, pe, ke, temp, press);
	}

	// Close the energy output file 
	fclose(e_file);

	return 0;
}
