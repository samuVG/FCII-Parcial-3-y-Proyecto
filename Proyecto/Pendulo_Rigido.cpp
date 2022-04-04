#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<vector>
#include "Pendulo_Rigido.h"

using namespace std;

//Runge--Kutta 4
//Resuelve un sistema con varias ecuaciones diferenciales 
vector<vector<double>> RK4(vector<double> (*f)(double, vector<double>, double, double, double, double, double), double a, vector<double> y0, double b, double n, double A, double B, double R1, double m1, double m2)
{
	int d = y0.size();
	
	//declaro 4 vectores
	vector<double> k1, k2, k3, k4;
	//lleno los vectores de d ceros
	for(int i=0; i<d; i++)
	{
		k1.push_back(0.);
		k2.push_back(0.);
		k3.push_back(0.);
		k4.push_back(0.);
	}
	
	//creo una matriz d+1xn, usando el contructor de la clase vector
	//contendra las posiciones y velocidades angulares de las masas 1 y 2
	//y la ultima fila de la matriz contiene los pasos de tiempo
	vector<vector<double>> Y(d+1, vector<double>(n));
	
	//paso de tiempo
	double h = (b-a)/n;
	
	//contador de ciclo while
	int ti=0;
	
	while(ti<n)
	{
		//se almacenan las coords al tiempo ti
		for(int i=0; i<d; i++){
			Y[i][ti] = y0[i];
		}
		
		//se almacena el tiempo
		Y[d][ti] = a;
		
		//declaro vectores auxiliares
		vector<double> aux1,aux2,aux3;
		//lleno los vectores con d ceros
		for(int i=0; i<d; i++)
		{
			aux1.push_back(0.);
			aux2.push_back(0.);
			aux3.push_back(0.);
		}
		
		// *** Runge-Kutta 4 metodo ***
		
		//se calculan los vectores k del metodo
		//calculo vector k1
		for(int i=0; i<d; i++)
		{	
			k1[i] = h*f(a,y0,A,B,R1,m1,m2)[i];
			aux1[i] = y0[i] + 0.5*k1[i];
		}
		
		//calculo vector k2
		for(int i=0; i<d; i++)
		{	
			k2[i] = h*f( a+0.5*h , aux1,A,B,R1,m1,m2 )[i];
			aux2[i] = y0[i] + 0.5*k2[i];
		}
		
		//calculo vector k3
		for(int i=0; i<d; i++)
		{	
			k3[i] = h*f( a+0.5*h , aux2,A,B,R1,m1,m2 )[i];
			aux3[i] = y0[i] + k3[i];
		}
		
		//calculo vector k4
		for(int i=0; i<d; i++)
		{	
			k4[i] = h*f( a+h , aux3,A,B,R1,m1,m2 )[i];
		}
		
		//actualizacion del vector inicial
		for(int i=0; i<d; i++)
		{	
			y0[i] = y0[i] + (k1[i] + 2*(k2[i]+k3[i]) +k4[i] )/6;
		}
		
		//actualizacion del tiempo y contador del while
		a += h;
		ti++;
	}
	
	//retorna matriz donde cada fila corresponde a la soluciones de la ED x1,x1',x2,x2',...,t; cada fila con n elementos
	return Y; 
}


void Pendulo_Doble_Rigido(vector<double> (*f)(double, vector<double>, double, double, double, double, double), double a, vector<double> y0, double b, double n, double A, double B, double R1, double m1, double m2) 
{
	//Distancia del pivote 1 al pivote 2
	double L1 = sqrt( pow(A/2,2) + pow(B,2) );
	
	//Angulos
	double beta = atan(A/(2*B));
	double eta = atan(A/B);
	double phi = 2*beta; //angulo entre L1 y R1
	
	//Distancia del pivote 2 al centro de masa del rectangulo 2
	double R2 = sqrt( pow(A/2,2)+pow(B/2,2) );
	
	//Distancia del centro de masa del rectangulo 1 al centro de masa X1
	double r = sqrt( pow(R1,2) + pow(B/2,2) - R1*B*cos(phi-beta) );
	
	//Momento de inercia de cada rectangulo
	double I1 = m1*( pow(A,2)+pow(B,2) )/12 + m1*pow(r,2);
	double I2 = m2*( pow(A,2)+pow(B,2) )/12; 
	
	//Soluciones de la ED x1,x1',x2,x2',...,t
	vector<vector<double>> Y = RK4(*(f), a, y0, b, n,A,B,R1,m1,m2); 
	
	//declaro los vectores correspondientes a los puntos moviles en la simulacion
	vector<double> x1, y1, x2, y2, p2x, p2y, v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y, v5x, v5y, v6x, v6y;
	//lleno los vectores con n ceros
	for(int i=0; i<n; i++)
	{
		x1.push_back(0);
		y1.push_back(0);
		x2.push_back(0);
		y2.push_back(0);
		p2x.push_back(0);
		p2y.push_back(0);
		v1x.push_back(0);
		v1y.push_back(0);
		v2x.push_back(0);
		v2y.push_back(0);
		v3x.push_back(0);
		v3y.push_back(0);
		v4x.push_back(0);
		v4y.push_back(0);
		v5x.push_back(0);
		v5y.push_back(0);
		v6x.push_back(0);
		v6y.push_back(0);
	}
	
	for(int i=0; i<n; i++)
	{
		//coords en x,y del Centro de masa de los rectangulos 1 y 2
		x1[i] = R1*sin(Y[0][i]); 
		y1[i] = -R1*cos(Y[0][i]);
		x2[i] = L1*sin(Y[0][i] + phi) + R2*sin(Y[2][i]);
		y2[i] = -L1*cos(Y[0][i] + phi) - R2*cos(Y[2][i]);
		
		//Posiciones de los vertices de los rectangulos y el pivote 2 en cada paso de tiempo
		v1x[i] = (A/2)*cos(Y[0][i]+phi-beta);
		v1y[i] = (A/2)*sin(Y[0][i]+phi-beta);
		
		v2x[i] = -v1x[i]; 
		v2y[i] = -v1y[i]; 
		
		v3x[i] = (B-(A/2)*(1/tan(Y[0][i]+phi-beta)))*sin(Y[0][i]+phi-beta);
		v3y[i] = ((A/2)*(1/tan(Y[0][i]+phi-beta))-B)*cos(Y[0][i]+phi-beta) - A/(2*sin(Y[0][i]+phi-beta));
		
		p2x[i] = L1*sin(Y[0][i]+phi);
		p2y[i] = -L1*cos(Y[0][i]+phi);
		
		v4x[i] = p2x[i] + A*cos(Y[2][i]-eta);
		v4y[i] = p2y[i] + A*sin(Y[2][i]-eta);
		
		v5x[i] = p2x[i] + 2*R2*sin(Y[2][i]);
		v5y[i] = p2y[i] - 2*R2*cos(Y[2][i]);
		
		v6x[i] = p2x[i] + B*sin(Y[2][i]-eta);
		v6y[i] = p2y[i] - B*cos(Y[2][i]-eta); 
	}
	
	//Creacion del archivo de texto con los valores para hacer la simulacion
	ofstream datos;
	
	datos.open("Datos_Simulacion.txt",ios::out); 
	
	if(datos.fail())
   {
		cout<<"No se pudo crear el archivo";
		exit(1);
  	}
  	
  	datos<<setw(15)<<left<<"x1"<<setw(15)<<"y1"<<setw(15)<<"x2"<<setw(15)<<"y2"<<setw(15)<<"v1x"
  		  <<setw(15)<<"v1y"<<setw(15)<<"v2x"<<setw(15)<<"v2y"<<setw(15)<<"v3x"<<setw(15)<<"v3y"
  		  <<setw(15)<<"v4x"<<setw(15)<<"v4y"<<setw(15)<<"v5x"<<setw(15)<<"v5y"<<setw(15)<<"v6x"
  		  <<setw(15)<<"v6y"<<setw(15)<<"p2x"<<setw(15)<<"p2y"<<setw(15)<<"Teta1"<<setw(15)<<"Teta2"
  		  <<setw(15)<<"Omega1"<<setw(15)<<"Omega2"<<endl;
  		  
	for (int i = 0; i < n; i++)
   {
   	datos<<setw(15)<<left<<x1[i]<<setw(15)<<y1[i]<<setw(15)<<x2[i]<<setw(15)<<y2[i]<<setw(15)<<v1x[i]
  		  <<setw(15)<<v1y[i]<<setw(15)<<v2x[i]<<setw(15)<<v2y[i]<<setw(15)<<v3x[i]<<setw(15)<<v3y[i]
  		  <<setw(15)<<v4x[i]<<setw(15)<<v4y[i]<<setw(15)<<v5x[i]<<setw(15)<<v5y[i]<<setw(15)<<v6x[i]
  		  <<setw(15)<<v6y[i]<<setw(15)<<p2x[i]<<setw(15)<<p2y[i]<<setw(15)<<Y[0][i]<<setw(15)<<Y[2][i]
  		  <<setw(15)<<Y[1][i]<<setw(15)<<Y[3][i]<<endl;
	}
   
   datos.close();
	
}











