#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<string>

#include "DiffFinNonLin.h"

using namespace std;

double Dfy(double (*f)(double, double, double), double x, double y, double dy)
{
	const double h=1.0e-5;
	
	return (f(x,y+h,dy) - f(x,y-h,dy)) / (2*h);	
}

double Dfdy(double (*f)(double, double, double), double x, double y, double dy)
{
	const double h=1.0e-5;
	
	return (f(x,y,dy+h) - f(x,y,dy-h)) / (2*h);
}

void Diferencias_Finitas(double (*f)(double,double,double), double (*y_exacta)(double), double a, double b, double alpha, double beta, int N, double TOL, int M)
{
	double w[N+2]; //arreglo con las soluciones de y
	double v[N+2]; //arreglo de la ecuacion Av = d 
   double x[N+2]; //arreglo con los valores de x
   double dy[N+2]; //arreglo con las soluciones de la primera derivada de y
   double u[N+2]; //arreglo U de la factorizacion de Crout A = LU
   double z[N+2]; //arreglo de la ecuacion Lz = d
   double l; //valor de la componente de la diagonal de L 
   double A; //coeficiente de la diagonal principal de la matriz tridiagonal
   double B; //coeficiente de la diagonal secundaria superior de la matriz tridiagonal
   double C; //coeficiente de la diagonal secundaria inferior de la matriz tridiagonal
   double d; //coeficiente del vector d de la ecuacion Av = d
   double h; //paso entre valores de x
  
	//PASO 1
   h = (b - a) / (N+1);
   
   x[0] = a;   
   x[N+1] = b;
   
   w[0] = alpha;
   w[N+1] = beta;
   
   z[0] = 0.; //se inicializa el primer elemento de estos arreglos
   u[0] = 0.;
   
   //PASO 2
	//primera aproximacion a la solucion   
   for (int i = 1; i < N+1; i++)
   {
   	w[i] = alpha + i*( (beta-alpha) /(b-a) )*h; 
   }
   
   //PASO 3
   int k = 1; //contador 
   
   //PASO 4
   //dy[0] = (w[1] - alpha)/(2.*h); 
   while (k < M + 1)
   {   
      //PASO 5,8    
		x[1] = a + h;
      dy[1] = (w[2] - alpha)/(2.*h);
      A = 2 + pow(h, 2) * Dfy( *(f), x[1], w[1], dy[1]);
      B = h/2 * Dfdy( *(f), x[1], w[1], dy[1]) - 1.;
      d = -( 2.*w[1] - w[2] - alpha + pow(h, 2)*f(x[1], w[1], dy[1]) );
      
      l = A;
      u[1] = B/l;
      z[1] = d/l;
         
      //PASO 6,9
      for (int i = 2; i < N; i++)
      {
      	x[i] = a + i*h;
         dy[i] = ( w[i+1] - w[i-1] ) / (2.*h);
         A = 2. + pow(h, 2) * Dfy( *(f), x[i], w[i], dy[i]);
         B = h/2. * Dfdy(*(f), x[i], w[i], dy[i]) - 1.;
         C = -1 - h/2. * Dfdy(*(f), x[i], w[i], dy[i]);
         d = -( 2.*w[i] - w[i+1] - w[i-1] + pow(h, 2)*f(x[i], w[i], dy[i]) );

         l = A - C*u[i-1];
         u[i] = B/l;
         z[i] = ( d - C*z[i-1] )/l;
		}
		
		//PASO 7,10
      x[N] = b - h;
      dy[N] = (beta - w[N-1]) / (2*h);
      A = 2. + pow(h, 2) * Dfy( *(f), x[N], w[N], dy[N]);
      C = -1 - h/2 * Dfdy(*(f), x[N], w[N], dy[N]);
      d = -( 2.*w[N] - w[N-1] - beta + pow(h, 2) * f(x[N], w[N], dy[N]) );
      
      l = A - C*u[N-1];
      z[N] = ( d - C*z[N-1] )/l;
			
      //PASO 11
      v[N] = z[N];
      w[N] += v[N];
        
      //PASO 12
      for (int i = N-1; i > 0 ; i--)
      {
      	v[i] = z[i] - u[i]*v[i+1];
         w[i] += v[i];
      }
        
      
      double magnitudv = 0.; 
      //para trabajar con los items originales de v y no modificarlos
      for (const auto &item : v ) 
      {
      	magnitudv += pow(item, 2);
      } 
		
		//PASO 13
      if (sqrt(magnitudv) <= TOL )
      {   
      	ofstream datos;
      	
      	string FinDiff="FinDiff_", Alpha_="a=", Beta_="_b=", txt=".txt";
      	string NombreArchivo = FinDiff + Alpha_ + to_string(a) + Beta_ + to_string(b)+txt;
      	
      	datos.open(NombreArchivo,ios::out); 
      	
      	if(datos.fail())
      	{
				cout<<"No se pudo crear el archivo";
				exit(1);
  			}
  			
  			cout<<endl;
  			cout<<setw(15)<<left<<"Xi"<<setw(15)<<"Wi"<<setw(15)<<"y_exacta"<<endl;
  			datos<<endl;
  			datos<<setw(15)<<left<<"Xi"<<setw(15)<<"Wi"<<setw(15)<<"y_exacta"<<endl;	  
			
			//PASO 14
         for (int i = 0; i < N+2; i++)
         {  
         	cout<<setw(15)<<left<<x[i]<<setw(15)<<w[i]<<setw(15)<<y_exacta(x[i])<<endl;
         	datos<<setw(15)<<left<<x[i]<<setw(15)<<w[i]<<setw(15)<<y_exacta(x[i])<<endl;
         }
         datos.close();
         //PASO 15
         break;      
		}
      //PASO 16
   	k++;    
	}
	//PASO 17
	cout<<"Se ha alcanzado el numero maximo de iteraciones. Si los resultados no son acertados pruebe aumentando este numero o evaluando los criterios de convergencia del metodo."<<endl;
}
   
   
   
   
   
   
   
 
