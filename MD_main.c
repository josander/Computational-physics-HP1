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

	// Initiation of variables 
	lattice_param = 4.05; // Units: [Å]
	timestep = 0.1; // 
	nbr_of_timesteps = 1000;
	nbr_of_atoms = 256;
	Nx = 4, Ny = 4, Nz = 4;
	m = 0.00279636665; // Metal units [ev/Å]

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
	for(i = 0; i < nbr_of_atoms; i++){
		for(j = 0; j < 3; j++){
			a[i][j] = f[i][j]/m;
		}
	}

	// Calculate of initial energies
	pe = get_energy_AL(q, Nx*lattice_param, nbr_of_atoms);
	ke = get_ke(v, nbr_of_atoms, m);
	energy = pe + ke;

	printf("E: %F \t Pe: %F \n", energy, pe);

	// Make a file to save the energies in
	FILE *e_file;
	e_file = fopen("energy.data","w");

	// Save the initial energies in the file
	fprintf(e_file,"%.5f \t %e \t %e \t %e \n", 0.0, energy, pe, ke);

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
		for(i = 0; i < nbr_of_atoms; i++){
			for(j = 0; j < 3; j++){
				a[i][j] = f[i][j]/m;
			}
		}

		// v(t+dt)
		for (j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
		    		v[j][n] += timestep * 0.5 * a[j][n];
			}
		} 

		// Calcutaion of the pe, ke and total energy
		pe = get_energy_AL(q, Nx*lattice_param, nbr_of_atoms);
		ke = get_ke(v, nbr_of_atoms, m);
		energy = pe + ke;
	
		// Print the average energy data to output file
		fprintf(e_file,"%.5f \t %e \t %e \t %e \n", i*timestep, energy, pe, ke);
	}

	// Close the energy output file 
	fclose(e_file);

	return 0;
}
