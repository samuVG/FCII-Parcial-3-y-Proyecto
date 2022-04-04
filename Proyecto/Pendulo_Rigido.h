#include<vector>

using namespace std;

//f:ecuaciones diferenciales acopladas, vector de condiciones iniciales, intervalo de integracion, numero de pasos de tiempo
vector<vector<double>> RK4(vector<double> (*f)(double, vector<double>, double, double, double, double, double), double, vector<double>, double, double, double, double, double, double, double);

//dimensiones rectangulos, distancia del pivote 1 al centro de masa, masa de los rectangulos
void Pendulo_Doble_Rigido(vector<double> (*f)(double, vector<double>, double, double, double, double, double), double, vector<double>, double, double, double, double, double, double, double);
