#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <complex>

#define N 200
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define n_ciclos 50
#define lambda 0.5
#define pasos 1000

using namespace std;

double calcular_norma(complex<double> phi[N+1]);
void generar(double& k, double& s, double V[N+1], complex<double> phi[N+1], complex<double> gamma[N], complex<double> alpha[N]);
void calcular_beta(double s, complex<double> gamma[N], complex<double> beta[N], complex<double> phi[N+1], complex<double> alpha[N], complex<double> xi[N+1]);

int main()
{

    double k, s, V[N+1];
    complex<double> phi[N+1], gamma[N], alpha[N], beta[N], xi[N+1];
    int j, l;
    ofstream datos;

    datos.open("schrodinger_data.dat");

    generar(k,s,V,phi,gamma,alpha);

    beta[N-1]=xi[N]=xi[0]=0;

    for(l=0;l<pasos;l++)
    {
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
    x0=N/4;
    sigma=N/16;

    norma=0;
    for (j=0;j<=N;j++)
    {
        if(j<(2*N/5) || j>(3*N/5))
        {
            V[j]=0;
        }
        else
        {
            V[j]=lambda*k*k;
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