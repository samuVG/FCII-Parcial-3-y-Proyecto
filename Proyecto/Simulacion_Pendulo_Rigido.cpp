#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<vector>
#include "Pendulo_Rigido.h"

using namespace std;

vector<double> f(double a, vector<double> y0, double A, double B, double R1, double m1, double m2)
{
	//Gravedad
	double g = 9.8;
	//vector donde almaceno los f1,f2,...
	vector<double> F;
	
	//Distancia del pivote 1 al pivote 2
	double L1 = sqrt( pow(A/2,2) + pow(B,2) );
	
	//Angulos
	double beta = atan(A/(2*B));
	double phi = 2*beta; //angulo entre L1 y R1
	
	//Distancia del pivote 2 al centro de masa del rectangulo 2
	double R2 = sqrt( pow(A/2,2)+pow(B/2,2) );
	
	//Distancia del centro de masa del rectangulo 1 al centro de masa X1
	double r = sqrt( pow(R1,2) + pow(B/2,2) - R1*B*cos(phi-beta) );
	
	//Momento de inercia de cada rectangulo
	double I1 = m1*( pow(A,2)+pow(B,2) )/12 + m1*pow(r,2);
	double I2 = m2*( pow(A,2)+pow(B,2) )/12;
	
	//Primera derivada de teta1, posicion angular para el centro de masa del rectangulo 1
	double f1 = y0[1];
	F.push_back(f1);
	
	//Segunda derivada de teta1, segunda derivada de la posicion angular para el centro de masa del rectangulo 1
	double num1 = -2*g*m1*R1*(I2+m2*pow(R2,2))*sin(y0[0]);
	double num2 = -L1*m2*g*(2*I2+m2*pow(R2,2))*sin(y0[0]+phi);
	double num3 = -L1*m2*R2*g*m2*R2*sin(y0[0]-2*y0[2]+phi);
	double num4 = -L1*m2*R2*2*(pow(y0[3],2)*(I2+m2*pow(R2,2))+pow(y0[1],2)*L1*m2*R2*cos(y0[0]-y0[2]+phi))*sin(y0[0]-y0[2]+phi);
	double den = 2*I2*L1*L1*m2+2*I2*m1*R1*R1+pow(L1*m2*R2,2)+2*m1*m2*pow(R1*R2,2)+2*I1*(I2+m2*R2*R2)-pow(L1*m2*R2,2)*cos(2*(y0[0]-y0[2]+phi));
	double f2 = (num1+num2+num3+num4)/den;
	F.push_back(f2);
	
	//Primera derivada de teta2, posicion angular para el centro de masa del rectangulo 2
	double f3 = y0[3];
	F.push_back(f3);
	
	//Segunda derivada de teta2, segunda derivada de la posicion angular para el centro de masa del rectangulo 2
	num1 = -m2*R2*g*(2*I1+L1*L1*m2+2*m1*R1*R1)*sin(y0[2]);
	num2 = m2*R2*L1*(g*m1*R1*sin(y0[2]-phi)+2*pow(y0[1],2)*(I1+L1*L1*m2+m1*R1*R1)*sin(y0[0]-y0[2]+phi));
	num3 = m2*R2*L1*(pow(y0[3],2)*L1*m2*R2*sin(2*(y0[0]-y0[2]+phi))+g*m1*R1*sin(2*y0[0]-y0[2]+phi)+g*L1*m2*sin(2*y0[0]-y0[2]+2*phi));
	den = 2*I2*L1*L1*m2+2*I2*m1*R1*R1+pow(L1*m2*R2,2)+2*m1*m2*pow(R1*R2,2)+2*I1*(I2+m2*R2*R2)-pow(L1*m2*R2,2)*cos(2*(y0[0]-y0[2]+phi));
	double f4 = (num1+num2+num3)/den;
	F.push_back(f4);
	
	return F;
}


int main()
{

	//Intervalo de integracion
	double a = 0, b = 80;
	
	//Numero de pasos de tiempo
	double n = 1000;	
	
	//Vector de condiciones iniciales del sistema
	double teta1 = M_PI/6, teta2 = 0.; //posiciones iniciales
	double w1 = 0., w2 = 0.; //velocidades angulares iniciales
	
	vector<double> y0 = {teta1, w1, teta2, w2};
	
	//Dimensiones de los rectangulos
	double A = 6, B = 12;
	
	//Distancia del pivote 1 al centro de masa del rectangulo 1
	double R1 = 7;
	
	//masas de los rectangulos
	double m1 = 1, m2 = 1;
	
	Pendulo_Doble_Rigido(*(f), a, y0, b, n, A, B, R1, m1, m2);
	
	return 0;
}






