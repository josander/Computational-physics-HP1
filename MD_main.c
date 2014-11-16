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


/* Main program */
int main()
{

    	/* Declaration of variables */
	double lattice_param;
	double timestep;
	double nbr_of_timesteps;
	double nbr_of_atoms;
	int Nx, Ny, Nz;
     	double random_value;
	int i, j;

	/* Initiation of variables */
	lattice_param = 4.05; // Units: [Å]
	timestep = 0.1;
	nbr_of_timesteps = 1000;
	nbr_of_atoms = 256;
	Nx = 4, Ny = 4, Nz = 4;

	/* Declaration of matrixes and arrays */
	double positions[4*Nx*Ny*Nz][3];
	double q[nbr_of_atoms], v[nbr_of_atoms], a[nbr_of_atoms];
	double energy[nbr_of_atoms];
	double pe[nbr_of_atoms];
	double ke[nbr_of_atoms];

	/* Initiation of the fcc lattice of Al-atoms */
	init_fcc(positions, Nx, lattice_param);

	// Introduction of displacements, about 5 % of the lattice spacing
     	srand(time(NULL)); // Should only be called once
     	random_value = (double) rand() / (double) RAND_MAX;

	/* Initiation of the velocities */
	for(i = 0; i < nbr_of_atoms; i++){
		v[i] = 0.0;
	}

	// Use the Al-potential to calculate the initial acceleration

	/* Calculation of initial energies */
	//calc_pe();
	//calc_ke();
	//energy = pe + ke;

	/* Make a file to save the energies in*/
	FILE *e_file;
	e_file = fopen("energy.data","w");

	/* Save the initial energies in the file*/
	fprintf(e_file,"%.5f \t %e \t %e \t %e \n", 0.0, energy, pe, ke);

	/* Time evolution according to the velocity Verlet algorithm */
	for (i = 1; i < nbr_of_timesteps + 1; i++){
		/* v(t+dt/2) */
		for (j = 0; j < nbr_of_atoms; j++){
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* q(t+dt) */
		for (j = 0; j < nbr_of_atoms; j++){
		    q[j] += timestep * v[j];
		}

		/* a(t+dt) */
		calc_acc(a, q, m, kappa, nbr_of_atoms);

		/* v(t+dt) */
		for (j = 0; j < nbr_of_atoms; j++){
		    v[j] += timestep * 0.5 * a[j];
		} 


		/* Calcutaion of the pe, ke and total energy */
		// calc_pe();
		// calc_ke();
		// energy = pe + ke;
	
		/* Print the average energy data to output file */
		fprintf(e_file,"%.5f \t %e \t %e \t %e \n", i*timestep, energy, pe, ke);
	}

	/* Close the energy output file */
	fclose(e_file);

	return 0;
}
