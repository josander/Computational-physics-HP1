/*
alpotential.c
 
Created by Anders on 2013-03-14.
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*Parameters for the AL EAM potential */
#define PAIR_POTENTIAL_ROWS 18 
const double pair_potential[90] = {2.0210, 2.2730, 2.4953, 2.7177, 2.9400, 3.1623, 3.3847, 3.6070, 3.8293, 4.0517, 4.2740, 4.4963, 4.7187, 4.9410, 5.1633, 5.3857, 5.6080, 6.0630, 2.0051, 0.7093, 0.2127, 0.0202, -0.0386, -0.0492, -0.0424, -0.0367, -0.0399, -0.0574, -0.0687, -0.0624, -0.0492, -0.0311, -0.0153, -0.0024, -0.0002, 0, -7.2241, -3.3383, -1.3713, -0.4753, -0.1171, 0.0069, 0.0374, 0.0122, -0.0524, -0.0818, -0.0090, 0.0499, 0.0735, 0.0788, 0.0686, 0.0339, -0.0012, 0, 9.3666, 6.0533, 2.7940, 1.2357, 0.3757, 0.1818, -0.0445, -0.0690, -0.2217, 0.0895, 0.2381, 0.0266, 0.0797, -0.0557, 0.0097, -0.1660, 0.0083, 0, -4.3827, -4.8865, -2.3363, -1.2893, -0.2907, -0.3393, -0.0367, -0.2290, 0.4667, 0.2227, -0.3170, 0.0796, -0.2031, 0.0980, -0.2634, 0.2612, -0.0102, 0};
    

#define ELECTRON_DENSITY_ROWS 15
const double electron_density[75] = {2.0210, 2.2730, 2.5055, 2.7380, 2.9705, 3.2030, 3.4355, 3.6680, 3.9005, 4.1330, 4.3655, 4.5980, 4.8305, 5.0630, 6.0630, 0.0824, 0.0918, 0.0883, 0.0775, 0.0647, 0.0512, 0.0392, 0.0291, 0.0186, 0.0082, 0.0044, 0.0034, 0.0027, 0.0025, 0.0000, 0.0707, 0.0071, -0.0344, -0.0533, -0.0578, -0.0560, -0.0465, -0.0428, -0.0486, -0.0318, -0.0069, -0.0035, -0.0016, -0.0008, 0, -0.1471, -0.1053, -0.0732, -0.0081, -0.0112, 0.0189, 0.0217, -0.0056, -0.0194, 0.0917, 0.0157, -0.0012, 0.0093, -0.0059, 0, 0.0554, 0.0460, 0.0932, -0.0044, 0.0432, 0.0040, -0.0392, -0.0198, 0.1593, -0.1089, -0.0242, 0.0150, -0.0218, 0.0042, 0};

#define EMBEDDING_ENERGY_ROWS 13
const double embedding_energy[65] = {0, 0.1000, 0.2000, 0.3000, 0.4000, 0.5000, 0.6000, 0.7000, 0.8000, 0.9000, 1.0000, 1.1000, 1.2000, 0, -1.1199, -1.4075, -1.7100, -1.9871, -2.2318, -2.4038, -2.5538, -2.6224, -2.6570, -2.6696, -2.6589, -2.6358, -18.4387, -5.3706, -2.3045, -3.1161, -2.6175, -2.0666, -1.6167, -1.1280, -0.4304, -0.2464, -0.0001, 0.1898, 0.2557, 86.5178, 44.1632, -13.5018, 5.3853, -0.3996, 5.9090, -1.4103, 6.2976, 0.6785, 1.1611, 1.3022, 0.5971, 0.0612, -141.1819, -192.2166, 62.9570, -19.2831, 21.0288, -24.3978, 25.6930, -18.7304, 1.6087, 0.4704, -2.3503, -1.7862, -1.7862};


/* Evaluates the spline in x. */

double splineEval(double x, const double *table,int m) {
	/* int m = mxGetM(spline), i, k;*/
	int i, k;
    
    /*double *table = mxGetPr(spline);*/
	double result;
    
	int k_lo = 0, k_hi = m;
    
	/* Find the index by bisection. */
	while (k_hi - k_lo > 1) {
		k = (k_hi + k_lo) >> 1;
		if (table[k] > x)
			k_hi = k;
		else
			k_lo = k;
	}
    
	/* Switch to local coord. */
	x -= table[k_lo];
    
	/* Horner's scheme */
	result = table[k_lo + 4*m];
	for (i = 3; i > 0; i--) {
		result *= x;
		result += table[k_lo + i*m];
	}
    
	return result;
}

/* Evaluates the derivative of the spline in x. */

double splineEvalDiff(double x, const double *table, int m) {
	/*int m = mxGetM(spline), i, k;
	double *table = mxGetPr(spline);
	*/
    int i, k;
    double result;
    
	int k_lo = 0, k_hi = m;
    
	/* Find the index by bisection. */
	while (k_hi - k_lo > 1) {
		k = (k_hi + k_lo) >> 1;
		if (table[k] > x)
			k_hi = k;
		else
			k_lo = k;
	}
    
	/* Switch to local coord. */
	x -= table[k_lo];
    
	/* Horner's scheme */
	result = 3*table[k_lo + 4*m];
	for (i = 3; i > 1; i--) {
		result *= x;
		result += (i-1)*table[k_lo + i*m];
	}
    
	return result;
}

/* Returns the forces */
void get_forces_AL(double forces[][3], double positions[][3], double cell_length, int nbr_atoms)
{
    int i, j;
    double cell_length_inv, cell_length_sq;
    double rcut, rcut_sq;
    double densityi, dens, drho_dr, force;
    double dUpair_dr;
    double sxi, syi, szi, sxij, syij, szij, rij,  rij_sq;
    
    double *sx = malloc(nbr_atoms * sizeof (double));
    double *sy = malloc(nbr_atoms * sizeof (double));
    double *sz = malloc(nbr_atoms * sizeof (double));
    double *fx = malloc(nbr_atoms * sizeof (double));
    double *fy = malloc(nbr_atoms * sizeof (double));
    double *fz = malloc(nbr_atoms * sizeof (double));
    
    double *density = malloc(nbr_atoms * sizeof (double));
    double *dUembed_drho = malloc(nbr_atoms * sizeof (double));
    
    rcut = 6.06;
    rcut_sq = rcut * rcut;
    
    cell_length_inv = 1 / cell_length;
    cell_length_sq = cell_length * cell_length;
    
    for (i = 0; i < nbr_atoms; i++){
        sx[i] = positions[i][0] * cell_length_inv;
        sy[i] = positions[i][1] * cell_length_inv;
        sz[i] = positions[i][2] * cell_length_inv;
    }
    
    for (i = 0; i < nbr_atoms; i++){
        density[i] = 0;
        fx[i] = 0;
        fy[i] = 0;
        fz[i] = 0;
    }
    
    for (i = 0; i < nbr_atoms; i++) {
        /* Periodically translate coords of current particle to positive quadrants */
		sxi = sx[i] - floor(sx[i]);
		syi = sy[i] - floor(sy[i]);
		szi = sz[i] - floor(sz[i]);
        
        densityi = density[i];
		
		/* Loop over other atoms. */
		for (j = i + 1; j < nbr_atoms; j++) {
            /* Periodically translate atom j to positive quadrants and calculate distance to it. */
			sxij = sxi - (sx[j] - floor(sx[j]));
			syij = syi - (sy[j] - floor(sy[j]));
			szij = szi - (sz[j] - floor(sz[j]));
            
            /* Periodic boundary conditions. */
			sxij = sxij - (int)floor(sxij + 0.5);
			syij = syij - (int)floor(syij + 0.5);
			szij = szij - (int)floor(szij + 0.5);
            
            /* squared distance between atom i and j */
			rij_sq = cell_length_sq * (sxij*sxij + syij*syij + szij*szij);
            
            /* Add force and energy contribution if distance between atoms smaller than rcut */
			if (rij_sq < rcut_sq) {
                rij = sqrt(rij_sq);
                dens = splineEval(rij, electron_density, ELECTRON_DENSITY_ROWS);
                densityi += dens;
                density[j] += dens;
            }
        }
        density[i] = densityi;
    }
    
    /* Loop over atoms to calculate derivative of embedding function
     and embedding function. */
	for (i = 0; i < nbr_atoms; i++) {
		dUembed_drho[i] = splineEvalDiff(density[i], embedding_energy, EMBEDDING_ENERGY_ROWS);
	}
    
    /* Compute forces on atoms. */
	/* Loop over atoms again :-(. */
    
    for (i = 0; i < nbr_atoms; i++) {
        /* Periodically translate coords of current particle to positive quadrants */
		sxi = sx[i] - floor(sx[i]);
		syi = sy[i] - floor(sy[i]);
		szi = sz[i] - floor(sz[i]);
        
        densityi = density[i];
		
		/* Loop over other atoms. */
		for (j = i + 1; j < nbr_atoms; j++) {
            /* Periodically translate atom j to positive quadrants and calculate distance to it. */
			sxij = sxi - (sx[j] - floor(sx[j]));
			syij = syi - (sy[j] - floor(sy[j]));
			szij = szi - (sz[j] - floor(sz[j]));
            
            /* Periodic boundary conditions. */
			sxij = sxij - (int)floor(sxij + 0.5);
			syij = syij - (int)floor(syij + 0.5);
			szij = szij - (int)floor(szij + 0.5);
            
            /* squared distance between atom i and j */
			rij_sq = cell_length_sq * (sxij*sxij + syij*syij + szij*szij);
            
            /* Add force and energy contribution if distance between atoms smaller than rcut */
			if (rij_sq < rcut_sq) {
                rij = sqrt(rij_sq);
                dUpair_dr = splineEvalDiff(rij, pair_potential, PAIR_POTENTIAL_ROWS);
				drho_dr = splineEvalDiff(rij, electron_density, ELECTRON_DENSITY_ROWS);
                
                /* Add force contribution from i-j interaction */
				force = -(dUpair_dr + (dUembed_drho[i] + dUembed_drho[j])*drho_dr) / rij;
				fx[i] += force * sxij * cell_length;
				fy[i] += force * syij * cell_length;
				fz[i] += force * szij * cell_length;
				fx[j] -= force * sxij * cell_length;
				fy[j] -= force * syij * cell_length;
				fz[j] -= force * szij * cell_length;
			}
        }
    }
    
    for (i = 0; i < nbr_atoms; i++){
        forces[i][0] = fx[i];
        forces[i][1] = fy[i];
        forces[i][2] = fz[i];
    }
    
    free(sx); free(sy); free(sz); sx = NULL; sy = NULL; sz = NULL;
    free(fx); free(fy); free(fz); fx = NULL; fy = NULL; fz = NULL;
    free(density); density = NULL;
    free(dUembed_drho); dUembed_drho = NULL;
    
}

/* Returns the potential energy */
double get_energy_AL(double positions[][3], double cell_length, int nbr_atoms)
{
    int i, j;
    double cell_length_inv, cell_length_sq;
    double rcut, rcut_sq;
    double energy;
    double densityi, dens;
    double sxi, syi, szi, sxij, syij, szij, rij,  rij_sq;
    
    double *sx = malloc(nbr_atoms * sizeof (double));
    double *sy = malloc(nbr_atoms * sizeof (double));
    double *sz = malloc(nbr_atoms * sizeof (double));
    
    double *density = malloc(nbr_atoms * sizeof (double));
    
    rcut = 6.06;
    rcut_sq = rcut * rcut;
    
    cell_length_inv = 1 / cell_length;
    cell_length_sq = cell_length * cell_length;
    
    for (i = 0; i < nbr_atoms; i++){
        sx[i] = positions[i][0] * cell_length_inv;
        sy[i] = positions[i][1] * cell_length_inv;
        sz[i] = positions[i][2] * cell_length_inv;
    }
    
    for (i = 0; i < nbr_atoms; i++){
        density[i] = 0;
    }
    
    energy = 0;
    
    for (i = 0; i < nbr_atoms; i++) {
        /* Periodically translate coords of current particle to positive quadrants */
		sxi = sx[i] - floor(sx[i]);
		syi = sy[i] - floor(sy[i]);
		szi = sz[i] - floor(sz[i]);
        
        densityi = density[i];
		
		/* Loop over other atoms. */
		for (j = i + 1; j < nbr_atoms; j++) {
            /* Periodically translate atom j to positive quadrants and calculate distance to it. */
			sxij = sxi - (sx[j] - floor(sx[j]));
			syij = syi - (sy[j] - floor(sy[j]));
			szij = szi - (sz[j] - floor(sz[j]));
            
            /* Periodic boundary conditions. */
			sxij = sxij - (int)floor(sxij + 0.5);
			syij = syij - (int)floor(syij + 0.5);
			szij = szij - (int)floor(szij + 0.5);
            
            /* squared distance between atom i and j */
			rij_sq = cell_length_sq * (sxij*sxij + syij*syij + szij*szij);
            
            /* Add force and energy contribution if distance between atoms smaller than rcut */
			if (rij_sq < rcut_sq) {
                rij = sqrt(rij_sq);
                dens = splineEval(rij, electron_density, ELECTRON_DENSITY_ROWS);
                densityi += dens;
                density[j] += dens;
                
                /* Add energy contribution from i-j interaction */
                energy += splineEval(rij, pair_potential, PAIR_POTENTIAL_ROWS);
                
            }
        }
        density[i] = densityi;
    }
    
    /* Loop over atoms to calculate derivative of embedding function
     and embedding function. */
	for (i = 0; i < nbr_atoms; i++) {
		energy += splineEval(density[i], embedding_energy, EMBEDDING_ENERGY_ROWS);
	}
    
    free(sx); free(sy); free(sz); sx = NULL; sy = NULL; sz = NULL;
    free(density); density = NULL;
    
    return(energy);

}

/* Returns the virial */ 
double get_virial_AL(double positions[][3], double cell_length, int nbr_atoms)
{
    int i, j;
    double cell_length_inv, cell_length_sq;
    double rcut, rcut_sq;
    double virial;
    double densityi, dens, drho_dr, force;
    double dUpair_dr;
    double sxi, syi, szi, sxij, syij, szij, rij, rij_sq;
    
    double *sx = malloc(nbr_atoms * sizeof (double));
    double *sy = malloc(nbr_atoms * sizeof (double));
    double *sz = malloc(nbr_atoms * sizeof (double));
    
    double *density = malloc(nbr_atoms * sizeof (double));
    double *dUembed_drho = malloc(nbr_atoms * sizeof (double));
    
    rcut = 6.06;
    rcut_sq = rcut * rcut;
    
    cell_length_inv = 1 / cell_length;
    cell_length_sq = cell_length * cell_length;
    
    for (i = 0; i < nbr_atoms; i++){
        sx[i] = positions[i][0] * cell_length_inv;
        sy[i] = positions[i][1] * cell_length_inv;
        sz[i] = positions[i][2] * cell_length_inv;
    }
    
    for (i = 0; i < nbr_atoms; i++){
        density[i] = 0;
    }
    
    for (i = 0; i < nbr_atoms; i++) {
        /* Periodically translate coords of current particle to positive quadrants */
		sxi = sx[i] - floor(sx[i]);
		syi = sy[i] - floor(sy[i]);
		szi = sz[i] - floor(sz[i]);
        
        densityi = density[i];
		
		/* Loop over other atoms. */
		for (j = i + 1; j < nbr_atoms; j++) {
            /* Periodically translate atom j to positive quadrants and calculate distance to it. */
			sxij = sxi - (sx[j] - floor(sx[j]));
			syij = syi - (sy[j] - floor(sy[j]));
			szij = szi - (sz[j] - floor(sz[j]));
            
            /* Periodic boundary conditions. */
			sxij = sxij - (int)floor(sxij + 0.5);
			syij = syij - (int)floor(syij + 0.5);
			szij = szij - (int)floor(szij + 0.5);
            
            /* squared distance between atom i and j */
			rij_sq = cell_length_sq * (sxij*sxij + syij*syij + szij*szij);
            
            /* Add force and energy contribution if distance between atoms smaller than rcut */
			if (rij_sq < rcut_sq) {
                rij = sqrt(rij_sq);
                dens = splineEval(rij, electron_density, ELECTRON_DENSITY_ROWS);
                densityi += dens;
                density[j] += dens;
            }
        }
        density[i] = densityi;
    }
    
    /* Loop over atoms to calculate derivative of embedding function
     and embedding function. */
	for (i = 0; i < nbr_atoms; i++) {
		dUembed_drho[i] = splineEvalDiff(density[i], embedding_energy, EMBEDDING_ENERGY_ROWS);
	}
    
    /* Compute forces on atoms. */
	/* Loop over atoms again :-(. */
    
    virial = 0;
    
    for (i = 0; i < nbr_atoms; i++) {
        /* Periodically translate coords of current particle to positive quadrants */
		sxi = sx[i] - floor(sx[i]);
		syi = sy[i] - floor(sy[i]);
		szi = sz[i] - floor(sz[i]);
        
        densityi = density[i];
		
		/* Loop over other atoms. */
		for (j = i + 1; j < nbr_atoms; j++) {
            /* Periodically translate atom j to positive quadrants and calculate distance to it. */
			sxij = sxi - (sx[j] - floor(sx[j]));
			syij = syi - (sy[j] - floor(sy[j]));
			szij = szi - (sz[j] - floor(sz[j]));
            
            /* Periodic boundary conditions. */
			sxij = sxij - (int)floor(sxij + 0.5);
			syij = syij - (int)floor(syij + 0.5);
			szij = szij - (int)floor(szij + 0.5);
            
            /* squared distance between atom i and j */
			rij_sq = cell_length_sq * (sxij*sxij + syij*syij + szij*szij);
            
            /* Add force and energy contribution if distance between atoms smaller than rcut */
			if (rij_sq < rcut_sq) {
                rij = sqrt(rij_sq);
                dUpair_dr = splineEvalDiff(rij, pair_potential, PAIR_POTENTIAL_ROWS);
				drho_dr = splineEvalDiff(rij, electron_density, ELECTRON_DENSITY_ROWS);
                
                /* Add virial contribution from i-j interaction */
				force = -(dUpair_dr + (dUembed_drho[i] + dUembed_drho[j])*drho_dr) / rij;
                
                virial += force * rij_sq;
			}
        }
    }
    
    virial /= 3.0;
    
    free(sx); free(sy); free(sz); sx = NULL; sy = NULL; sz = NULL;
    free(density); density = NULL;
    free(dUembed_drho); dUembed_drho = NULL;
    
    return(virial);
    
}