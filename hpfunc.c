/*
hpfunc.c
Contains function for homeproblem 1/b

By Jossan och Svensson

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Randomly displaces the atoms in position[][3] by 5% of the lattice parameter
void rand_disp(double position[][3] ,double lattice_param ,int N)
{
	int i, j;
	double sum;	
	double ran[3];
	double length_corr;	
	for(i = 0; i < N; i++){
		sum = 0;		
		for(j = 0; j < 3; j++){
			ran[j] = rand() - 0.5 ; // ran[j] goes from -0.5 to 0.5
			sum =+ ran[j]*ran[j];
		}
		sum = sqrt(sum);
		length_corr = (N*0.05)/sum;
		
		for(j = 0; j < 3; j++){
			position[i][j] =+ ran[j] * length_corr;
		}

	}  	
		  
}
