#include "DiffFinNonLin.h"
#include<cmath>

using namespace std;

double f(double x, double y, double dy)
{		
	return - (pow(y,5)+2*dy/x);
}

double theta_exacta(double x)
{
	return 1/sqrt(1+pow(x,2)/3);
}

int main()
{
	double a = 0;
   double b = 5;
   double alpha = 1;
   double beta = sqrt(21)/14;
   double N = 100, M = 1e6; //h=0.1
   double TOL = 1e-8; 
    
   Diferencias_Finitas(*(f), *(theta_exacta), a, b, alpha, beta, N, TOL, M);
   
	return 0;
}
