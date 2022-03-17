#include "DiffFinNonLin.h"
#include<cmath>

using namespace std;

double f1(double x, double y, double dy)
{
    return (32 + 2 * pow(x, 3) - y * dy) / 8;
}

double y_exacta1(double x)
{
	return pow(x,2) + 16/x;
}

double f2(double x, double y, double dy)
{
    return pow(y, 3) - y * dy;
}

double y_exacta2(double x)
{
	return pow(x+1,-1);
}

int main()
{
	double a1 = 1., a2 = 1., a3 = 0.;
   double b1 = 3., b2 = 2., b3 = 5.;
   double alpha1 = 17.0, alpha2 = 1./2., alpha3 = 1;
   double beta1 = 43.0/3.0, beta2 = 1./3., beta3 = sqrt(21)/14;
   int N = 19, N3 = 100; //h=0.1
   int M = 10, M3 = 1e6;
   double TOL = 1e-8; 
    
   Diferencias_Finitas(*(f1), *(y_exacta1), a1, b1, alpha1, beta1, N, TOL, M);
   Diferencias_Finitas(*(f2), *(y_exacta2), a2, b2, alpha2, beta2, N, TOL, M);
   
	return 0;
}








