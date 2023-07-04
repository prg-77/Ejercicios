#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <complex>

#define N 400
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define n_ciclos 50
#define lambda 0.5
#define pasos 2000
#define experimentos 1
#define n_max 20 //Numero maximo de barreras
#define n_barreras 20 //Numero de barreras

using namespace std;

double calcular_norma(complex<double> phi[N+1]);
void generar(double& k, double& s, double V[N+1], complex<double> phi[N+1], complex<double> gamma[N], complex<double> alpha[N]);
void calcular_beta(double s, complex<double> gamma[N], complex<double> beta[N], complex<double> phi[N+1], complex<double> alpha[N], complex<double> xi[N+1]);
double calcular_PD(complex<double> phi[N+1]);
int calcular_maximo(double PD[pasos]);
double calcular_x(complex<double> phi[N+1]);

int main()
{

    double k, s, V[N+1], PD[pasos], p, merr;
    complex<double> phi[N+1], gamma[N], alpha[N], beta[N], xi[N+1];
    int j, l, m, mT, mi[experimentos];
    ofstream datos;

    datos.open("schrodinger_data.dat");

    mT=0;

    srand(time(NULL));

    for(m=0;m<experimentos;m++)
    {
        generar(k,s,V,phi,gamma,alpha);

        beta[N-1]=xi[N]=xi[0]=0;

        for(l=0;l<pasos;l++)
        {
            
            PD[l]=calcular_PD(phi);
            
            if(l%3==0)
            {
                for(j=0;j<=N;j++)
                {
                    datos << j << ", " << abs(phi[j])*abs(phi[j]) << ", " << V[j] << ", " << real(phi[j])*real(phi[j]) << ", " << imag(phi[j])*imag(phi[j]) << ", " <<  calcular_norma(phi) <<  endl;
                }
                datos << endl;
            }
            
            calcular_beta(s,gamma,beta,phi,alpha,xi);

        }

        cout << "Maximo local de PD: " << calcular_maximo(PD) << endl;
        cout << "Valor de PD en el maximo local: " << PD[calcular_maximo(PD)] << endl;

        p=1.0*rand()/RAND_MAX;

        if(p<PD[calcular_maximo(PD)])
        {
            mT=mT+1;
            mi[m]=1;
        }
        else
        {
            mi[m]=0;
        }
    }

    merr=0;
    for(m=0;m<experimentos;m++)
    {
        merr=merr+(mi[m]-1.0*mT/experimentos)*(mi[m]-1.0*mT/experimentos);
    }
    merr=sqrt(merr/experimentos)/sqrt(experimentos);  // Media en un intervalo de confianza del 68%

    cout << "Probabilidad de transmision: " << 1.0*mT/experimentos << " +/- " << merr << endl;

    datos.close();

    return 0;   
}


double calcular_norma(complex<double> phi[N+1])
{
    int j;
    double norma;
    norma=0;
    for(j=0;j<=N;j++)
    {
        norma=norma+abs(phi[j])*abs(phi[j]);
    }
    return sqrt(norma);
}

void generar(double& k, double& s, double V[N+1], complex<double> phi[N+1], complex<double> gamma[N], complex<double> alpha[N])
{
    double norma, x0, sigma;
    int j,l;

    k=2.0*PI*n_ciclos/N;
    s=1.0/(4*k*k);
    x0=N/100;
    sigma=N/200;

    norma=0;
    for (j=0;j<=N;j++)
    {
        V[j]=0;

        for(l=0;l<n_barreras;l++)
        {
            if(j>=(N/(2*n_max+1)+(2*l*N/(2*n_max+1))) & j<=((2*N/(2*n_max+1))+(2*l*N/(2*n_max+1))))
            {
                V[j]=lambda*k*k;
            }
        }
    
        phi[j]=pow(exp(1.0i),k*j)*exp(-(j-x0)*(j-x0)/(2*sigma*sigma));
        norma=norma+abs(phi[j])*abs(phi[j]);
    }
    phi[0]=phi[N]=0;

    for(j=0;j<=N;j++)
    {
        phi[j]=phi[j]/sqrt(norma);
    }

    alpha[N-1]=0;
    for(j=N-1;j>0;j--)
    {
        gamma[j]=1.0/(-2-V[j]+2.0i/s+alpha[j]);
        alpha[j-1]=-1.0*gamma[j];
    }

    return;
}

void calcular_beta(double s, complex<double> gamma[N], complex<double> beta[N], complex<double> phi[N+1], complex<double> alpha[N], complex<double> xi[N+1])
{
    int j;
    
    for(j=N-1;j>0;j--)
    {
        beta[j-1]=gamma[j]*(4.0i*phi[j]/s-beta[j]);
    }
    for(j=0;j<(N-1);j++)
    {
        xi[j+1]=alpha[j]*xi[j]+beta[j];
        phi[j+1]=xi[j+1]-phi[j+1];
    }

    return;

}

double calcular_PD(complex<double> phi[N+1])
{
    int j;
    double PD;
    PD=0;
    for(j=(N-(N/(2*n_max+1)));j<=N;j++)  // No se si habria que poner mas amplio el intervalo
    {
        PD=PD+abs(phi[j])*abs(phi[j]);
    }
    return PD;
}

// Para calcular el primer maximo local de PD se recorre el array PD hasta que el siguiente valor sea menor que el anterior
int calcular_maximo(double PD[pasos])
{
    int j;
    bool mayor;

    mayor=true;
    j=0;
    
    while(mayor==true)
    {
        if(PD[j+1]>=PD[j])
        {
            j++;
        }
        else
        {
            mayor=false;
        }
    }
    return j;
}

// Cálculo del observable x
double calcular_x(complex<double> phi[N+1])
{
    int j;
    double x;
    x=0;
    for(j=0;j<=N;j++)
    {
        x=x+j*abs(phi[j])*abs(phi[j]);
    }
    return x;
}
