#include <cstdio>
#include <ctime>
#include "PSOC.hpp"

const unsigned int NUM_PARTICULA = 75;
const unsigned char NUM_PPP = 7;
const float C1 = 2.0;
const float C2 = 2.0;
const float InerMax = 0.9;
const float InerMin = 0.4;
const double LimInferiorR = 1000;
const double LimSuperiorR = 1000000;
const double LimInferiorC = 0.000000001;
const double LimSuperiorC = 0.000001;
const float VmaxR = 1.0;
const float VminR = -1.0;
//const float VmaxC = 0.0000000002;
//const float VminC = -0.0000000002;
const float VmaxC = 0.000000000001;
const float VminC = -0.000000000001;
const unsigned int NUM_MAX_ITER = 10000; 

using namespace std;

int main()
{
	unsigned int t = 0, n, m; 
	float inercia = 0.0;
	float error;	
	//Numero de Particula, Numero de Parametros, Xmin, Xmax, Numero maximo de iteraciones
	//PSOC Circuit;
	PSOC Circuit(NUM_PARTICULA , NUM_PPP, C1, C2, LimInferiorR, LimSuperiorR, LimSuperiorC, LimInferiorC);
	srand((unsigned int)time(NULL));
	Circuit.InicializarEnjambre();
	//Circuit.MuestraEnjambre();
	Circuit.EvaluacionInicial();
	//Circuit.MuestraEnjambre();
	Circuit.ActualizarMejores();
	printf("\n\n La mejor es: \n");
	Circuit.MuestraMejor();
	//error = 1;
	//Sintonizacion del Algoritmo
	//Explosion del Enjambre
	//while((t<NUM_MAX_ITER)&&(error>0.000001))
	//scanf("%d", &m);
	while(t<NUM_MAX_ITER)
	{
		inercia = InerMax - (((InerMax - InerMin)/NUM_MAX_ITER) * t);
		//printf("%f\n", inercia);
		//scanf("%d", &n);
		Circuit.ActualizarVelocidad(inercia, VmaxR, VminR, VmaxC, VminC);
		Circuit.ActualizarPosicion();
		Circuit.EvaluarEnjambre();
		Circuit.ActualizarMejores();
		//printf("\n\n Primera Particula: \n");
		//Test.MuestraEnjambre();
		printf("\n\n Iteracion: %i\n", t);		
		//Test.MuestraEnjambre();
		printf("\n\n %i Mejor Particula: \n",t);
		Circuit.MuestraMejor();
		//error = 50 - Test.GetBestFit();
		t++;
	}
	return 0;
}