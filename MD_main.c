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
	int i, j;

	// Initiation of variables 
	lattice_param = 4.05; // Units: [Å]
	timestep = 0.1;
	nbr_of_timesteps = 1000;
	nbr_of_atoms = 256;
	Nx = 4, Ny = 4, Nz = 4;

	// Declaration of matrixes and arrays 
	double q[4*Nx*Ny*Nz][3];
	double v[nbr_of_atoms], a[nbr_of_atoms];
	double energy[nbr_of_atoms];
	double pe[nbr_of_atoms];
	double ke[nbr_of_atoms];

	// Initiation of the fcc lattice of Al-atoms 
	init_fcc(q, Nx, lattice_param);

	// Introduction of displacements, about 5 % of the lattice spacing
     	srand(time(NULL)); // Should only be called once
     	//random_value = (double) rand() / (double) RAND_MAX;

	
	rand_disp(q, lattice_param, nbr_of_atoms);
	/*
	for(i = 0; i < 10; i++){
		printf("%F \t %F \t %F \n", q[i][0], q[i][1], q[i][2]);
	}
	*/
	// Initiation of the velocities 
	for(i = 0; i < nbr_of_atoms; i++){
		v[i] = 0.0;
	}

	// Use the Al-potential to calculate the initial acceleration

	/* Calculation of initial energies */
	//get_energy_AL(q, Nx*lattice_param, nbr_of_atoms);
	//calc_ke();
	//energy = pe + ke;

	/* Make a file to save the energies in*/
	FILE *e_file;
	e_file = fopen("energy.data","w");

	/* Save the initial energies in the file*/
	fprintf(e_file,"%.5f \t %e \t %e \t %e \n", 0.0, energy, pe, ke);



	/* Close the energy output file */
	fclose(e_file);

	return 0;
}
