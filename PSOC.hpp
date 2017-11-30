#include<cstdio>
#include<cstdlib>
#include<cmath>

using namespace std;

typedef struct{
	double *Xi;
	double *Vi;
	double *Pi;
	double rXFit;
	double rPFit;
	double SwR1;
	double SwR2;
	double SqR1;
	double SqR2;
	double SqR4;
	double SqC35;
	double Eg;
	double Ewp;
	double Eq;
	double rPSwR1;
	double rPSwR2;
	double rPSqR1;
	double rPSqR2;
	double rPSqR4;
	double rPSqC35;
	double rPEg;
	double rPEwp;
	double rPEq;
}PARTICULA;

void FuncionObjetivo(PARTICULA *Pi, double Gfi, double Wpfi, double Qpfi);
double FuncionFitness(PARTICULA *Pi, double EmGi, double Emwpi, double EmQi, double maxEgGenera, double maxEwpGenera, double maxEqGenera);

class PSOC
{
private:
	//unsigned int NumMaxIteraciones;
	unsigned int NumParticulas;
	unsigned char NumParametros;
	float C1;
	float C2;
	double XminR;
	double XmaxR;
	double XminC;
	double XmaxC;
	double EmG;
	double Emwp;
	double EmQ;
	double Gf;
	double Wpf;
	double Qpf;
	PARTICULA *Enj;
	unsigned int IdMejorP;

public:
	PSOC(unsigned int Np, unsigned char dim, float Ci1, float Ci2, double XiniR, double XfinR, double XiniC, double XfinC);
	void InicializarEnjambre(void);
	void MuestraEnjambre(void);
	void EvaluarEnjambre(void);
	void EvaluacionInicial(void);
	void ActualizarMejores(void);
	void MuestraMejor(void);
	void ActualizarVelocidad(double Iner, float VmaxR, float VminR, float VmaxC, float VminC);
	void ActualizarPosicion(void);
	float GetBestFit();
	~PSOC();
};

double FuncionFitness(PARTICULA *Pi, double EmGi, double Emwpi, double EmQi, double maxEgGenera, double maxEwpGenera, double maxEqGenera)
{
	double r, aux1, aux2, aux3;
	double wj = (1/6);

	Pi->SwR1 = abs(Pi->SwR1);
	Pi->SwR2 = abs(Pi->SwR2);
	Pi->SqR1 = abs(Pi->SqR1);
	Pi->SqR2 = abs(Pi->SqR2);
	Pi->SqR4 = abs(Pi->SqR4);
	Pi->SqC35 = abs(Pi->SqC35);
	aux1 = abs(EmGi - Pi->Eg) / maxEgGenera;
	aux2 = abs(Emwpi - Pi->Ewp) / maxEwpGenera;
	aux3 = abs(EmQi - Pi->Eq) / maxEqGenera;
	r = (wj*Pi->SwR1) + (wj*Pi->SwR2) + (wj*Pi->SqR1) + (wj*Pi->SqR2) + (wj*Pi->SqR4) + (wj*Pi->SqC35);
	r += aux1 + aux2 + aux3;
	return r;
}

void FuncionObjetivo(PARTICULA *Pi, double Gfi, double Wpfi, double Qpfi)
{
	int n;
	double Ga, G, QP, R1, R2,Ra, R4, Rb, C3, C5, WP;
	
	//Se normaliza para evitar desbordamiento en la memoria.
	R1 = Pi->Xi[0]/1000;
	R2 = Pi->Xi[1]/1000;
	Ra = Pi->Xi[2]/1000;
	R4 = Pi->Xi[3]/1000;
	Rb = Pi->Xi[4]/1000;
	C3 = Pi->Xi[5]*1000;
	C5 = Pi->Xi[6]*1000;

	//Ga = 1 + (Pi->Xi[2]/Pi->Xi[4]); 
	Ga = 1+(Ra/Rb);
	// G = (Ga/(Pi->Xi[0]*Pi->Xi[6])) * (1/((1/(Pi->Xi[0]*Pi->Xi[6])) + (1/(Pi->Xi[3]*Pi->Xi[6])) + 
	// 	(1/(Pi->Xi[3]*Pi->Xi[5])) + ((1-Ga)/(Pi->Xi[1]*Pi->Xi[6]))));
	G = (Ga/(R1*C5))*(1/((1/(R1*C5))+(1/(R4*C5))+(1/(R4*C3))+((1-Ga)/(R2*C5))));
	//WP = sqrt( (  (Pi->Xi[0]+Pi->Xi[1]) / (Pi->Xi[0]*Pi->Xi[6]*Pi->Xi[3]*Pi->Xi[5]*Pi->Xi[1]) ) );
	WP = sqrt((R1+R2)/(R1*R2*R4*C5*C3));
	//QP = WP * (1/((1/(Pi->Xi[0]*Pi->Xi[6]))+(1/(Pi->Xi[3]*Pi->Xi[6]))+(1/(Pi->Xi[3]*Pi->Xi[5])) +
	//	((1-Ga)/(Pi->Xi[1]*Pi->Xi[6]))));
	QP = WP * (1 / ( (1 / (R1*C5) ) + (1 / (R4*C5) ) + (1 / (R4*C3) ) + ((1-Ga) / (R2*C5) ) ));
	//Pi->SwR1 = (-1*(Pi->Xi[1]))/(2*(Pi->Xi[0] + Pi->Xi[1]));
	Pi->SwR1 = (-1*(R2))/(2*(R1+R2));
	//Pi->SwR2 = (-1*(Pi->Xi[0]))/(2*(Pi->Xi[0] + Pi->Xi[1]));
	Pi->SwR2 = (-1*(R1))/(2*(R1+R2));
	//Pi->SqR1 = (Pi->Xi[1] * ((Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[3]) - (Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[1]) - 
	//	(Pi->Xi[6]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[5]*Pi->Xi[1]*Pi->Xi[3]) + (Pi->Xi[5]*Ga*Pi->Xi[0]*Pi->Xi[3]))) /
	//	(2 * (Pi->Xi[0]+Pi->Xi[1]) * ((Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[3]) +
	//	(Pi->Xi[6]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[5]*Pi->Xi[1]*Pi->Xi[3]) - (Pi->Xi[5]*Ga*Pi->Xi[0]*Pi->Xi[3])));
	Pi->SqR1 = (R2*((C3*R1*R4)-(C3*R1*R2)-(C5*R1*R2)+(C3*R2*R4)+(C3*Ga*R1*R4)))/(2*(R1+R2)*((C3*R1*R2)+(C3*R1*R4)+(C5*R1*R2)+(C3*R2*R4)-(C3*Ga*R1*R4)));
	//Pi->SqR2 = -1 * (((Pi->Xi[0]*(Pi->Xi[0]*Pi->Xi[1]*(Pi->Xi[5] + Pi->Xi[6]))) + (Pi->Xi[0]*Pi->Xi[3]*Pi->Xi[5]*(Ga - 1)) + 
	//	(Pi->Xi[1]*Pi->Xi[3]*Pi->Xi[5]*((2*Ga)-1))) / 
	//	(2 * (Pi->Xi[0]+Pi->Xi[1]) * ((Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[3]) +
	//	(Pi->Xi[6]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[5]*Pi->Xi[1]*Pi->Xi[3]) - (Pi->Xi[5]*Ga*Pi->Xi[0]*Pi->Xi[3]))));
	Pi->SqR2 = -1*((R1*((R1*R2*(C3+C5))+(R1*R4*C3*(Ga-1))+(R2*R4*C3*((2*Ga)-1))))/(2*(R1+R2)*((C3*R1*R2)+(C3*R1*R4)+(C5*R1*R2)+(C3*R2*R4)-(C3*Ga*R1*R4))));
	//Pi->SqR4 = (((Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[6]*Pi->Xi[0]*Pi->Xi[1])) / ((Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[1]) + 
	//	(Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[3]) + (Pi->Xi[6]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[5]*Pi->Xi[1]*Pi->Xi[3]) - 
	//	(Pi->Xi[5]*Ga*Pi->Xi[0]*Pi->Xi[3]))) - (1/2);
	Pi->SqR4 = ( ( (C3*R1*R2) + (C5*R1*R2) ) / ( (C3*R1*R2) + (C3*R1*R4) + (C5*R1*R2) + (C3*R2*R4) - (C3*Ga*R1*R4)) ) - (1/2);
	//Pi->SqC35 = ((Pi->Xi[6]*Pi->Xi[0]*Pi->Xi[1]) / ((Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[1]) + 
	//	(Pi->Xi[5]*Pi->Xi[0]*Pi->Xi[3]) + (Pi->Xi[6]*Pi->Xi[0]*Pi->Xi[1]) + (Pi->Xi[5]*Pi->Xi[1]*Pi->Xi[3]) + 
	//	(Pi->Xi[5]*Ga*Pi->Xi[0]*Pi->Xi[3]))) - (1/2);
	Pi->SqC35 = ( (C5*R1*R2) / ( (C3*R1*R2) + (C3*R1*R4) + (C5*R1*R2) + (C3*R2*R4) + (C3*Ga*R1*R4)) ) - (1/2);
	//
	//Calculo del Error Obtenido de los componentes seleccionados.
	Pi->Eg = abs((G - Gfi) / Gfi);
	Pi->Ewp = abs((WP - Wpfi) / Wpfi);
	Pi->Eq = abs((QP - Qpfi)/ Qpfi);
}

void PSOC::ActualizarPosicion(void)
{
	for(unsigned int i=0; i < NumParticulas; i++)
		for(unsigned int d=0; d < NumParametros; d++)
			Enj[i].Xi[d] += Enj[i].Vi[d];
}

void PSOC::ActualizarVelocidad(double Iner, float VmaxR, float VminR, float VmaxC, float VminC)
{
	float Y1, Y2;
	int n;

	for(unsigned int i=0; i < NumParticulas; i++)
		for(unsigned int d=0; d < NumParametros; d++)
			{
				Y1 = ((float)rand()/(float)RAND_MAX);
				Y2 = ((float)rand()/(float)RAND_MAX);
				//Enj[i].Vi[d] = (Enj[i].Vi[d]*Iner) + C1*Y1*(Enj[i].Pi[d] - Enj[i].Xi[d]) + C2*Y2*(Enj[IdMejorP].Pi[d] - Enj[i].Xi[d]);
				Enj[i].Vi[d] *= Iner + C1*Y1*(Enj[i].Pi[d] - Enj[i].Xi[d]) + C2*Y2*(Enj[IdMejorP].Pi[d] - Enj[i].Xi[d]);
				if(d < 5)
				{
					if(Enj[i].Vi[d]>VmaxR)
					 	Enj[i].Vi[d]=VmaxR;
					if(Enj[i].Vi[d]<VminR)
					 	Enj[i].Vi[d]=VminR;
					if(Enj[i].Vi[d]==0)
						Enj[i].Vi[d]=VmaxR;
				}
				else
				{
					//scanf("%d", &n);
					if(Enj[i].Vi[d]>VmaxC)
						Enj[i].Vi[d]=VmaxC;
				 	if(Enj[i].Vi[d]<-VmaxC)
				 	 	Enj[i].Vi[d]=-VmaxC;
				 	if(Enj[i].Vi[d]==0)
						Enj[i].Vi[d]=VmaxC;
				}
				//printf("%2.21f\n", Enj[i].Vi[d]);
				//if(Enj[i].Vi[d] == 0)
				//	scanf("%d", &n);
			}
	//scanf("%d", &n);		
}

void PSOC::ActualizarMejores(void)
{
	//Busca saber si la posicion en Xi es mejor que la que ya hay en Pi
	for(unsigned int i=0; i < NumParticulas; i++)
	{
		//Actualiza la mejor posicion de cada particula
		if (Enj[i].rXFit < Enj[i].rPFit)
			for(int k=0; k < NumParametros; k++)
				{
					Enj[i].Pi[k] = Enj[i].Xi[k];
					if(i == k)
					{
						Enj[i].rPFit = Enj[i].rXFit;
						Enj[i].rPSwR1 = Enj[i].SwR1;
						Enj[i].rPSwR2 = Enj[i].SwR2;
						Enj[i].rPSqR1 = Enj[i].SqR1;
						Enj[i].rPSqR2 = Enj[i].SqR2;
						Enj[i].rPSqR4 = Enj[i].SqR4;
						Enj[i].rPSqC35 = Enj[i].SqC35;
						Enj[i].rPEg = Enj[i].Eg;
						Enj[i].rPEwp = Enj[i].Ewp;
						Enj[i].rPEq = Enj[i].Eq;
					}
				}
		//Actualizar la mejor particula del Enjambre
		if (Enj[i].rXFit < Enj[IdMejorP].rPFit)
			IdMejorP = i;
	}
}

void PSOC::EvaluarEnjambre(void)
{
	double maxEg = 0.0, maxEwp = 0.0, maxEq = 0.0;
	//Para cada particula, evalua la funcion objetivo
	for(unsigned int i=0; i < NumParticulas; i++)
	{
		FuncionObjetivo(&Enj[i], Gf, Wpf, Qpf);
	}
	for(unsigned int j=0; j < NumParticulas-1; j++)
	{
		if(Enj[j].Eg >= Enj[j+1].Eg)
			maxEg = Enj[j].Eg;
		else
			maxEg = Enj[j+1].Eg;
		//Ciclo para obtener el error maximo en Ewp
		if(Enj[j].Ewp >= Enj[j+1].Ewp)
			maxEwp = Enj[j].Ewp;
		else
			maxEwp = Enj[j+1].Ewp;
		//Ciclo para obtener el error maximo en Eq
		if(Enj[j].Eq >= Enj[j+1].Eq)
			maxEq = Enj[j].Eq;
		else
			maxEq = Enj[j+1].Eq;
	}
	for(unsigned int k=0; k < NumParticulas; k++)
	{
		Enj[k].rXFit = FuncionFitness(&Enj[k], EmG, Emwp, EmQ, maxEg, maxEwp, maxEq);
	}
}

void PSOC::EvaluacionInicial(void)
{
	double maxEg = 0, maxEwp = 0, maxEq = 0;
	//Para cada particula, evalua la funcion objetivo
	for(unsigned int i=0; i < NumParticulas; i++)
	{
		FuncionObjetivo(&Enj[i], Gf, Wpf, Qpf);
	}
	for(unsigned int j=0; j < NumParticulas-1; j++)
	{
		if(Enj[j].Eg >= Enj[j+1].Eg)
			maxEg = Enj[j].Eg;
		else
			maxEg = Enj[j+1].Eg;
		//Ciclo para obtener el error maximo en Ewp
		if(Enj[j].Ewp >= Enj[j+1].Ewp)
			maxEwp = Enj[j].Ewp;
		else
			maxEwp = Enj[j+1].Ewp;
		//Ciclo para obtener el error maximo en Eq
		if(Enj[j].Eq >= Enj[j+1].Eq)
			maxEq = Enj[j].Eq;
		else
			maxEq = Enj[j+1].Eq;
	}
	//Para cada particula, evalua la funcion objetivo
	for(unsigned int k=0; k < NumParticulas; k++)
	{
		Enj[k].rXFit = FuncionFitness(&Enj[k], EmG, Emwp, EmQ, maxEg, maxEwp, maxEq);
		Enj[k].rPFit = Enj[k].rXFit;
		Enj[k].rPSwR1 = Enj[k].SwR1;
		Enj[k].rPSwR2 = Enj[k].SwR2;
		Enj[k].rPSqR1 = Enj[k].SqR1;
		Enj[k].rPSqR2 = Enj[k].SqR2;
		Enj[k].rPSqR4 = Enj[k].SqR4;
		Enj[k].rPSqC35 = Enj[k].SqC35;
		Enj[k].rPEg = Enj[k].Eg;
		Enj[k].rPEwp = Enj[k].Ewp;
		Enj[k].rPEq = Enj[k].Eq;
	}
}

void PSOC::MuestraMejor(void)
{
	//Muestra cada Particula en Pantalla
	printf("\nX%i: ", IdMejorP);
        for(int k=0; k < NumParametros; k++)
		printf("%+2.9f ", Enj[IdMejorP].Xi[k]);
	printf("\nV%i: ", IdMejorP);
	for(int k=0; k < NumParametros; k++)
		printf("%+2.12f ", Enj[IdMejorP].Vi[k]);
	printf("\nP%i: ", IdMejorP);
	for(int k=0; k < NumParametros; k++)
		printf("%+2.9f ", Enj[IdMejorP].Pi[k]);
	
	//Muestra los Xfitnes y los Pfitnes de la funcion
	printf("\nXFit%i: ", IdMejorP);
	printf("%+2.9f ", Enj[IdMejorP].rXFit);
	printf("\nPFit%i: ", IdMejorP);
	printf("%+2.9f ", Enj[IdMejorP].rPFit);

	printf("\nSwR1:");
	printf("%2.4f ", Enj[IdMejorP].rPSwR1);
	printf("\nSwR2:");
	printf("%2.4f ", Enj[IdMejorP].rPSwR2);
	printf("\nSqR1:");
	printf("%2.4f ", Enj[IdMejorP].rPSqR1);
	printf("\nSqR2:");
	printf("%2.4f ", Enj[IdMejorP].rPSqR2);
	printf("\nSqR4:");
	printf("%2.4f ", Enj[IdMejorP].rPSqR4);
	printf("\nSqC35:");
	printf("%2.4f ", Enj[IdMejorP].rPSwR1);
	printf("\nEg:");
	printf("%2.4f ", Enj[IdMejorP].rPEg);
	printf("\nEwp:");
	printf("%2.4f ", Enj[IdMejorP].rPEwp);
	printf("\nEq:\n");
	printf("%2.4f ", Enj[IdMejorP].rPEq);
}

void PSOC::MuestraEnjambre(void)
{
	int n;
	for(unsigned int i=0; i < NumParticulas; i++)
        {
        //Muestra cada Particula en Pantalla
		printf("\nX%i: ", i);
                for(int k=0; k < NumParametros; k++)
			printf("%+2.9f ", Enj[i].Xi[k]);
		printf("\nV%i: ", i);
		for(int k=0; k < NumParametros; k++)
			printf("%+2.9f ", Enj[i].Vi[k]);
		printf("\nP%i: ", i);
		for(int k=0; k < NumParametros; k++)
			printf("%+2.9f ", Enj[i].Pi[k]);
		
		//Muestra los Xfitnes y los Pfitnes de la funcion
		printf("\nXFit%i: ", i);
		printf("%+2.9f ", Enj[i].rXFit);
		printf("\nPFit%i: ", i);
		printf("%+2.9f ", Enj[i].rPFit);
        }
    //scanf("%d", &n);
}

void PSOC::InicializarEnjambre(void)
{
	double aux;
	for(int i=0; i < NumParticulas; i++)
    {
		//Inicializar la posicion de la particula Xi
		for(int k=0; k < NumParametros; k++)
		{
			//Calcula valor aletorio
			aux = (double)rand()/(double)RAND_MAX;
			//Calcula un valor aleatorio entre Rmin y Rmax
			if (k < 5)
				aux = XminR + aux*(XmaxR - XminR);
			else
				aux = XminC + aux*(XmaxC - XminC);
            Enj[i].Xi[k] = aux;
            Enj[i].Vi[k] = 0;
            Enj[i].Pi[k] = aux;
		}
	Enj[i].rXFit = 0.0;
	Enj[i].rPFit = 0.0;
	Enj[i].SwR1 = 0.0;
	Enj[i].SwR2 = 0.0;
	Enj[i].SqR1 = 0.0;
	Enj[i].SqR2 = 0.0;
	Enj[i].SqR4 = 0.0;
	Enj[i].SqC35 = 0.0;
	Enj[i].Eg = 0.0;
	Enj[i].Ewp = 0.0;
	Enj[i].Eq = 0.0;
	Enj[i].rPSwR1 = 0.0;
	Enj[i].rPSwR2 = 0.0;
	Enj[i].rPSqR1 = 0.0;
	Enj[i].rPSqR2 = 0.0;
	Enj[i].rPSqR4 = 0.0;
	Enj[i].rPSqC35 = 0.0;
	Enj[i].rPEg = 0.0;
	Enj[i].rPEwp = 0.0;
	Enj[i].rPEq = 0.0;
    }
}

PSOC::PSOC(unsigned int Np, unsigned char dim, float Ci1, float Ci2, double XiniR, double XfinR, double XiniC, double XfinC)
{
	//Define el numero de particulas del enjambre
	NumParticulas = Np;
	//Define el tamaño de los vectores de la particula
	NumParametros = dim;
	//Asigna las constantes C1 y C2
	C1 = Ci1;
	C2 = Ci2;
	//Asigna el rango del vector posicion en Resistencias 
	XminR = XiniR;
	XmaxR = XfinR;
	//Asigna el rango del vector posicion en Capacitores
	XminC = XiniC;
	XmaxC = XfinC;
	//Inicializa el Error maximo permitido de los parámetros funcionales del filtro
	EmG = 0.005;
	Emwp = 0.005;
	EmQ = 0.005;
	//
	Gf = 1;
	//Wpf = 4995.13*(0.159155/1);
	Wpf = 4995.13;
	Qpf = 1.11;
	//Inicializa el Indice de la mejor particula
	IdMejorP = 0;
	//Reserva memoria dinamica para el Enjambre
	Enj = new PARTICULA[NumParticulas];
	//Reserva memoria dinamica para los vectores Xi, Vi, Pi
	for(int i=0; i < NumParticulas; i++)
	{
		Enj[i].Xi = new double[NumParametros];
		Enj[i].Vi = new double[NumParametros];
		Enj[i].Pi = new double[NumParametros];
	}
}

PSOC::~PSOC()
{
	//Liberar memoria de los vectores
	for(int i=0; i < NumParticulas; i++)
        {
               delete [] Enj[i].Xi;
               delete [] Enj[i].Vi;
               delete [] Enj[i].Pi;
        }
	//Liberar memoria del enjambre
	delete [] Enj;
}