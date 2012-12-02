#include <stdio.h>
#include <stdlib.h>
#include "shared.h"

int progonka (double** koef_array , double* psi , double* f , const int n , const int n_0)
{
//	double *p = malloc(n * size of );
	double p[n+2+n_0];
	double q[n+2+n_0];
	
	if ( koef_array[n_0][n_0] != 0 )
		{
		p[n_0+1] = koef_array[n_0+1][n_0+2]/koef_array[n_0+1][n_0+1];
		q[n_0+1] = f[n_0]/koef_array[n_0+1][n_0+1];	
		}
	else 
		{
		printf("b0 = 0\n");
		exit (-2);
		}

int i;
	for (i = n_0 ; i <= n_0+n ; i++ )
	{	
		if ( koef_array[i][i] != (koef_array[i][i -1] *(-p[i-1]) )  )
		{
		p[i] = koef_array[i][i+1]/ (koef_array[i][i] - koef_array[i][i -1] *(0-p[i-1] ));
		q[i] = (koef_array[i][i-1] * p[i-1] - f[i-1])/ (koef_array[i][i] - koef_array[i][i -1] *( 0 -p[i-1] ));
		}
		else 
		{
			printf ("nesovmestnaya sistema\n");
			exit(-2);
		}
	}	

	
	psi[n_0+n] = (koef_array[n+n_0][n+n_0-1]*q[n_0+n] - f[n+n_0])/(koef_array[n+n_0][n+n_0] -
		koef_array[n+n_0][n+n_0-1]*p[n_0+n]); 
	for (i = n+n_0; i > n_0; i--)
	{
		psi[i-1] =  psi[i] * p[i] + q[i];

	}

	return 0;
}

int cycle_progonka (double** koef_array , double* psi , double* f , const int n , const int n_0)
{
double p[n_0+n-1];
double q[n_0+n-1];
double An[n_0+n-1];
int i;
//double matrix[n_0+n-1][n_0+n-1];
progonka (koef_array, p, f, n-1, n_0);
for (i = n_0; i <= n+n_0-2; i++)
	{
	An[i] = koef_array[i][n_0+n];
	
	} 
progonka (koef_array, q, An, n-1, n_0);
double sum = 0;
double sum2 = 0;
for (i = n_0; i< n+n_0-1; i++)
	{
	sum = sum + koef_array[n_0+n][i]*q[i];
	sum2 = sum2 + koef_array[n_0+n][i]*p[i];
	}
double g_fnAn[n+n_0-1];
for (i=n_0; i<= n_0+n-2; i++)
	{
	double  fn = (sum - f[n+n_0])/(sum2 - koef_array[n+n_0][n+n_0]);
	g_fnAn[i] = f[i] -(fn) * koef_array[i][n_0+n];
	}
progonka (koef_array, psi, g_fnAn, n-1, n_0);

	

return 0;
}

int psi_calc (struct field** array_b,struct field** space_array, double** koef_array, double** koef_array_half)
{
/*

	
	double f_array_half[hor_size][vert_size];
	double f_array[hor_size][vert_size];

//	double koef_array_half[hor_size][vert_size];	will define it seperately in function init_koef_array, and must be taken from it
//	double koef_array[hor_size][vert_size];
//	will accept it in a function

	double psi_line_half[vert_size];
	double psi_line[vert_size];
	int i,j;




	for (j = 1; j < vert_size; j++)
	{  
		for (i = 1; i < hor_size; i++)
		{ // initialysing f array for a half of a step
			f_array_half[i][j] = (-2) * alfa_x*array_b[i][j].psi -  ( array_b[i-1][j].psi - 2*array_b[i][j].psi
					+array_b[i+1][j].psi	 + array_b[i][j].omega*h_x*h_x) ;
		}
		// solving system of equesions
		cycle_progonka (koef_array_half, psi_line_half, (double**)f_array_half, n, n_0);
		//now we have psi_line_half initialized
		//  putting psi_line-half to arrya_b in the half of a step. it overwrites existing data form previous step
		//  because we don't use it anywhere and don't need it anymore 
		// so we don't need psi_matrix
		for  (i = 1; i < hor_size; i ++)
			{
				array_b[i][j].psi = psi_line_half[i];
			}	//psi_line -> psi_matrix (to  array_b ?)

	}
	



	for (i = 1; i < vert_size; i++) // including vert size??
	{  
		for (j = 1; j < hor_size ; j++)
		{ // initialysing f array for next part of step
			f_array[i][j] =  (-2)*alfa_z*array_b[i][j].psi - ( array_b[i][j-1].psi - 2*array_b[i][j].psi
					+array_b[i][j+1].psi  - h_z *0.5* sigma*(array_b[i][j+1].psi-array_b[i][j-1].psi));
				// from psi_matrix
			
		}
		//some edge conditions I've found!!!
		f_array[i][hor_size -1] = f_array[i][hor_size -1] - (1 + sigma*h_z)PSI;		
			
		//solving system of equations
		progonka (koef_array , psi_line  , (double**)f_array, n , n_0)

		//   putting psi to space_array in the end of a step

		for  (j = 1; j < hor_size; j++)
		{
			space_array[i][j].psi = psi_line_half[i];
				//psi_line -> psi_matrix (to  space_array ?)
		}
	}*/
return 0;
//copy psi_matrix to sp
}


int init_systems (void)
{
/*	double koef_array_half[hor_size][vert_size];
	double koef_array[hor_size][vert_size];
int i ;
	//definig koef_array_half
	koef_array_half[1][vert_size] = 1;
	koef_array_half[hor_size][1] = 1;
	koef_array_half[vert_size - 1][vert_size - 1] = -beta_x;
	for ( i = 1; i < vert_size - 1; i++ ) 
	{	
		koef_array_half[i][i+1] = 1;
		koef_array_half[i+1][i] = 1;
		koef_array_half[i][i] = -beta_x;
	}

//definig koef_array
	koef_array_half[vert_size - 1][vert_size - 1] = -beta_z;
	for ( i = 1; i < vert_size - 1; i++ ) 
	{	
		koef_array_half[i][i+1] = 1 + sigma*h_z/2;
		koef_array_half[i+1][i] = 1 - sigma*h_z/2;
		koef_array_half[i][i] = -beta_z;
	}

*/
	return 0;
}

int calc_omega_and_teta (struct field **space_array, struct field **array_b)
{
//	int i, j;

	/* double K_0;
	double h_x;
	double h_z;
	double sigma;
	double alfa[vert_size+3];
	double beta;
	double dt;   */
	
//	cp_array (array_b,array_b1);
//	cp_array (space_array,array_b);

/*
	for (i = 2; i < (vert_size - 1); i++){
		for (j = 2; j < (hor_size - 1); j++)
		{
			array_b[i][j].J11 = (1/(4*h_x*h_z))*( (array_b[i+1][j].psi - array_b[i-1][j].psi)*
					(array_b[i][j+1].omega - array_b[i][j-1].omega) - 
					(array_b[i][j+1].psi - array_b[i][j-1].psi)*
					(array_b[i+1][j].omega - array_b[i-1][j].omega)) + 
				sigma*(array_b[i+1][j].psi*array_b[i+1][j].omega - 
						array_b[i-1][j].psi*array_b[i-1][j].omega)/(2*h_x);

			array_b[i][j].J21 = (1/(4*h_x*h_z))*( array_b[i+1][j].psi*
					(array_b[i+1][j+1].omega - array_b[i+1][j-1].omega) - 
					array_b[i-1][j].psi*
					(array_b[i-1][j+1].omega - array_b[i-1][j-1].omega) - 
					array_b[i][j+1].psi*
					(array_b[i+1][j+1].omega - array_b[i-1][j+1].omega) + 
					array_b[i][j-1].psi*
					(array_b[i+1][j-1].omega + array_b[i-1][j-1].omega)) + 
				sigma*(array_b[i+1][j].psi*(array_b[i+1][j-1].omega + 
							array_b[i+1][j+1].omega ) - 
						array_b[i-1][j].psi*
						(array_b[i-1][j-1].omega +
						 array_b[i-1][j+1].omega))/(4*h_x);

			array_b[i][j].J12 = (1/(4*h_x*h_z))*( 0 - array_b[i+1][j].omega*
					(array_b[i+1][j+1].psi - array_b[i+1][j-1].psi) + 
					array_b[i-1][j].omega*
					(array_b[i-1][j+1].psi - array_b[i-1][j-1].psi) + 
					array_b[i][j+1].omega*
					(array_b[i+1][j+1].psi - array_b[i-1][j+1].psi) - 
					array_b[i][j-1].omega*
					(array_b[i+1][j-1].psi + array_b[i-1][j-1].psi)) +
				sigma*(array_b[i+1][j].omega*(array_b[i+1][j-1].psi + 
							array_b[i+1][j+1].psi ) - 
						array_b[i-1][j].omega*
						(array_b[i-1][j-1].psi + 
						 array_b[i-1][j+1].psi))/(4*h_x);

			array_b[i][j].J_omega = (array_b[i][j].J11 + array_b[i][j].J12 + 
					array_b[i][j].J21)/3;

			array_b[i][j].J_teta = (1/(2*h_x*h_z))*(( array_b[i -1][j].teta - 
					array_b[i][j].teta)*(array_b[i-1][j].psi - 
					array_b[i-1][j-1].psi) +
					(array_b[i][j].teta - array_b[i+1][j].teta)*
					(array_b[i][j].psi - array_b[i][j-1].psi) + 
					(array_b[i][j+1]. teta - array_b[i][j].teta)*
					(array_b[i][j].psi - array_b[i-1][j].psi) + 
					(array_b[i][j].teta - 
					 array_b[i][j-1].teta)*(array_b[i][j-1].psi - 
					 array_b[i-1][j-1].psi )+
					(sigma/(2*h_x))*( (array_b[i][j].teta - 
					array_b[i-1][j].teta)*
					(array_b[i-1][j-1].psi - array_b[i-1][j].psi) +	
					( array_b[i + 1][j].teta -array_b[i][j].teta)*
					(array_b[i][j-1].psi - 
					 array_b[i][j].psi)) );

			array_b[i][j].T_omega = K_0 * ( (array_b[i-1][j].omega - 2*array_b[i][j].omega 
						+ array_b[i+1][j].omega)/h_x/h_x  + 
					array_b[i][j-1].omega - 
					2*array_b[i][j].omega +array_b[i][j+1].omega) /h_z/h_z;
			array_b[i][j].T_teta =  K_0 * ( array_b[i-1][j].teta - 2*array_b[i][j].teta + 
					array_b[i+1][j].teta) /h_x/h_x + array_b[i][j -1].teta - 
				2*array_b[i][j].teta + array_b[i][j +1].teta /h_z/h_z
				+ (alfa[j] -alfa[j-1]) / h_z;

			array_b[i][j].F_omega = array_b[i][j].J_omega + beta/2/h_x * 
				( array_b[i+1][j].teta - 
				  array_b[i][j].teta + array_b[i+1][j+1].teta - 
				  array_b[i][j].teta)
				+  array_b[i][j].T_omega;
			array_b[i][j].F_teta = array_b[i][j].J_teta - alfa[j]/2/h_x * 
				(array_b[i][j].psi -array_b[i-1][j].psi+ array_b[i][j-1].psi - 
				 array_b[i-1][j-1].psi + array_b[i-1][j-1].psi) + array_b[i][j].T_teta;

			//NEW fields
			
			space_array[i][j].omega = array_b[i][j].omega + 
				1.5*dt *array_b[i][j].F_omega - 0.5 * dt * array_b1[i][j].F_omega;
			space_array[i][j].teta = array_b[i][j].teta + 
				1.5*dt *array_b[i][j].F_teta - 0.5 * dt * array_b1[i][j].F_teta; 


		}
	}; */
	return 0;	

}
