#include <stdio.h>
#include <stdlib.h>

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



int main (void)
{
double koef_array[10][10];
double f[10];
double psi[10];
progonka((double**)koef_array, psi, f, 10, 0);
return 0;


}
