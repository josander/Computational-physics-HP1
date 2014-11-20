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
#define nbr_of_atoms 256

// Main program 
int main()
{

    	// Declaration of variables 
	double lattice_param, cell_size;
	double timestep;
	int nbr_of_timesteps;
	int Nx, Ny, Nz;
     	double random_value;
	int i, j, k, n;
	double m;
	double energy, pe, ke;
	double tau_T, tau_P;
	double temp_eq;
	double press_eq;
	double kappa_P;
	int startCut;
	double msd;
	double vel;
	double self_diffusion;
	double meanF;
	int nbr_of_terms;


	// Initiation of variables 
	lattice_param = 4.05; // Units: [Å]
	timestep = 0.01; // [ps]
	nbr_of_timesteps = 2000;
	Nx = 4, Ny = 4, Nz = 4;
	m = 0.00279636665; // Metal units [ev/Å]
	temp_eq = 900 + 273.15; // Degree Celsius 
	press_eq = 6.324209 * pow(10, -7); // 1 Atm in eV/Å^3
	tau_T = timestep*100;
	tau_P = timestep*100;
	kappa_P = 2.21901454; //3.85 * pow(10, 9);/ // Liquid Aluminum Units: Å^3/eV
	cell_size = lattice_param*Nx;
	startCut = 1000;
	self_diffusion = 0;
	meanF = 0;

	// Declaration of matrixes and arrays 
	double q[4*Nx*Ny*Nz][3], v[nbr_of_atoms][3], a[nbr_of_atoms][3];
	double f[4*Nx*Ny*Nz][3];
	double *temp = malloc((nbr_of_timesteps+1) * sizeof(double));
	double *press = malloc((nbr_of_timesteps+1) * sizeof(double));
	double *corr_func_T = malloc((nbr_of_timesteps-startCut+1) * sizeof(double));
	double *corr_func_P = malloc((nbr_of_timesteps-startCut+1) * sizeof(double));
	double *MSD = malloc((nbr_of_timesteps+1) * sizeof(double));
	double *vel_corr_func = malloc((nbr_of_timesteps+1) * sizeof(double));
	double *spectral_func = malloc((nbr_of_timesteps+1) * sizeof(double));
	double *omega = malloc((nbr_of_timesteps+1) * sizeof(double));

	// Declaration of the Q-array and the V-array
	double *allElements = malloc((nbr_of_timesteps+1)*(nbr_of_atoms)*(3)*sizeof(double));
	double ***Q = malloc((nbr_of_timesteps+1)*sizeof(double **));
	double ***V = malloc((nbr_of_timesteps+1)*sizeof(double **));
	for(i = 0; i < nbr_of_timesteps +1; i++){
		Q[i] = malloc(nbr_of_atoms * sizeof(double *));
		V[i] = malloc(nbr_of_atoms * sizeof(double *));
		for(j = 0; j < nbr_of_atoms; j++){
			Q[i][j] = allElements + (i * nbr_of_atoms * 3) + (j * 3);
			V[i][j] = allElements + (i * nbr_of_atoms * 3) + (j * 3);
		}
	}

	// Initiation of Q and V
	for(i = 0; i < nbr_of_timesteps + 1; i++){
		for(j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
				Q[i][j][n] = 0.0;
				V[i][j][n] = 0.0;
			}
		}
	}

	// Initiation of corr_func
	for(i = 0; i <nbr_of_timesteps - startCut + 1; i++){
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

	// Save the inital displacement in Q
	for(j = 0; j < nbr_of_atoms; j++){
		for(n = 0; n < 3; n++){
			Q[0][j][n] = q[j][n];
			V[0][j][n] = v[j][n];
		}
	}

	// Get the initial Einstein frequency for the calculation of the spectral function
	for(j = 0; j < nbr_of_atoms; j++){
		meanF += (f[j][0]*f[j][0] + f[j][1]*f[j][1] + f[j][2]*f[j][2])/nbr_of_atoms;
	}
	omega[0] = sqrt(meanF/(3*m*K_B*temp[0]));

	// Make a file to save the energies and displacement of a particle in
	FILE *e_file;
	e_file = fopen("energy.data","w");

	FILE *d_file;
	d_file = fopen("displacement.data","w");

	// Save the initial energies in the file
	fprintf(e_file,"%.5f \t %e \t %e \t %e \t %F \t %e \n", 0.0, energy, pe, ke, temp[0], press[0]);
	fprintf(d_file,"%.5f \t %e \t %e \n", q[100][0], q[100][1], q[100][2]);

	// Time evolution according to the velocity Verlet algorithm (eqlib)
	for (i = 1; i < startCut ; i++){
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

		// Get the Einstein frequency for the calculation of the spectral function
		meanF = 0;
		for(j = 0; j < nbr_of_atoms; j++){
			meanF += (f[j][0]*f[j][0] + f[j][1]*f[j][1]+f[j][2]*f[j][2])/nbr_of_atoms;
		}
		omega[i] = sqrt(meanF/(3*m*K_B*temp[i]));

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

		// Save the displacement in Q
		for(j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
				Q[i][j][n] = q[j][n];
				V[i][j][n] = v[j][n];
			}
		}

		// Print the average energy data to output file
		fprintf(e_file,"%.5f \t %e \t %e \t %e \t %F \t %e \n", i*timestep, energy, pe, ke, temp[i], press[i]);
		fprintf(d_file,"%.5f \t %e \t %e \n", q[100][0], q[100][1], q[100][2]);
	}

	for (i = startCut; i < nbr_of_timesteps + 1; i++){ // No eqlib
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

		// Get the Einstein frequency for the calculation of the spectral function
		meanF = 0;
		for(j = 0; j < nbr_of_atoms; j++){
			meanF += (f[j][0]*f[j][0] + f[j][1]*f[j][1]+f[j][2]*f[j][2])/nbr_of_atoms;
		}
		omega[i] = sqrt(meanF/(3*m*K_B*temp[i]));

		
		// Calculate the pressure
		cell_size = lattice_param * Nx;
		press[i] = get_P(q, cell_size, nbr_of_atoms, temp[i]);
		
		// Calcutaion of the pe, ke and total energy
		pe = get_energy_AL(q, cell_size, nbr_of_atoms);
		pe = sqrt(pe*pe);
		ke = get_ke(v, nbr_of_atoms, m);
		energy = sqrt(pe*pe) + sqrt(ke*ke);

		// Print every 1000 timestep in the terminal
		if(i%1000 == 0){
			printf("%i av %i steg \n", i, nbr_of_timesteps);
		}

		// Save the displacement in Q
		for(j = 0; j < nbr_of_atoms; j++){
			for(n = 0; n < 3; n++){
				Q[i][j][n] = q[j][n];
				V[i][j][n] = v[j][n];
			}
		}

		// Print the average energy data to output file
		fprintf(e_file,"%.5f \t %e \t %e \t %e \t %F \t %e \n", i*timestep, energy, pe, ke, temp[i], press[i]);
		fprintf(d_file,"%.5f \t %e \t %e \n", q[100][0], q[100][1], q[100][2]);
	}


	// Get the correlation functions
	printf("TEMPERATURE: \n");
	get_corr_func(temp, corr_func_T, nbr_of_timesteps+1, startCut);
	printf("PRESSURE: \n");
	get_corr_func(press, corr_func_P, nbr_of_timesteps+1, startCut);

	// Write corr_func to data file
	FILE *c_file;
	c_file = fopen("correlation.data","w");

	for(i = 0; i < (nbr_of_timesteps-startCut+1); i++){
		fprintf(c_file,"%.5f \t %e \n", corr_func_T[i], corr_func_P[i]);
	}

	// Calculate the mean squared displacement and the velocity correlation function
	nbr_of_terms = 200 + 1;
	for(i = 0; i < nbr_of_terms; i++){
		for(j = 0; j < nbr_of_terms; j++){
			for(k = 0; k < nbr_of_atoms; k++){
				msd = sqrt((Q[i+j][k][0] - Q[i][k][0])*(Q[i+j][k][0] - Q[i][k][0]) + (Q[i+j][k][1] - Q[i][k][1])*(Q[i+j][k][1] - Q[i][k][1]) + (Q[i+j][k][2] - Q[i][k][2])*(Q[i+j][k][2] - Q[i][k][2]));
				MSD[j] += msd*msd/(nbr_of_terms*nbr_of_atoms);

				vel = (V[i+j][k][0] * V[i][k][0]) + (V[i+j][k][1] * V[i][k][1])+ (V[i+j][k][2] * V[i][k][2]);
				vel_corr_func[j] += vel/(nbr_of_terms*nbr_of_atoms);
			}
		}
	}

	// Calculate the spectral function
	get_spectral_func(vel_corr_func, omega, spectral_func, nbr_of_terms, timestep);

	// Calculate the self diffusion coefficient from the MSD. Valid method if q>>l and t>>tau
	self_diffusion = MSD[nbr_of_timesteps]/(6 * nbr_of_timesteps * timestep);

	// Print the self diffusion coefficient from the MSD in the terminal
	printf("Self diffusion coefficient from MSD: %e \n",self_diffusion);

	// Print the self diffusion coefficient from the vel-corr-func in the terminal
	printf("Self diffusion coefficient from vel-corr-func: %e \n",self_diffusion);

	// New file to print the MSD
	FILE *m_file;
	m_file = fopen("MSD.data","w");

	// Save the MSD-data
	for(j = 0; j < 200; j++){
		fprintf(m_file,"%.5f \t %e \t %e \t %e \n", timestep*j, MSD[j], vel_corr_func[j], spectral_func[j]);
	}

	// Close the energy output file 
	fclose(e_file);
	fclose(c_file);
	fclose(d_file);
	fclose(m_file);

	// Free allocated memory
	free(temp); free(press); free(corr_func_T); free(corr_func_P); free(Q); free(V); free(vel_corr_func); free(spectral_func); free(MSD); free(omega);
	temp = NULL; press = NULL; corr_func_T = NULL; corr_func_P = NULL; Q = NULL; V = NULL; vel_corr_func = NULL; spectral_func = NULL; MSD = NULL; omega = NULL; 

	return 0;
}
