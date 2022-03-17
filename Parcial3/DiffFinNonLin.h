#include<iostream>

using namespace std;


double Dfy(double (*f)(double, double, double), double, double, double);

double Dfdy(double (*f)(double, double, double), double, double, double);

void Diferencias_Finitas(double (*f)(double,double,double), double (*y_exacta)(double), double , double , double , double , int , double , int );
