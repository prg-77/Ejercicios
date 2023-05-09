#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <complex>

#define N 300
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define n_ciclos 75
#define lambda 0.8
#define pasos 200

using namespace std;

int main()
{

    double k,s, V[N+1], norma;
    complex<double> phi[N+1][pasos], gamma[N], alpha[N], beta[N], xi[N+1];
    int j,l;
    ofstream datos;

    k=2*PI*n_ciclos/N;
    s=1/(4*k*k);

    datos.open("schrodinger_data.dat");

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
        phi[j][0]=pow(exp(1.0i),k*j)*exp(-8.0*(4*j-N)*(4*j-N)/(N*N));
    }
    phi[0][0]=phi[N][0]=0;
    norma=0;
    for(j=0;j<=N;j++)
    {
        norma=norma+abs(phi[j][0])*abs(phi[j][0]);
    }
    for(j=0;j<=N;j++)
    {
        phi[j][0]=phi[j][0]/sqrt(norma);
    }
    

    alpha[N-1]=0;
    for(j=N-1;j>0;j--)
    {
        gamma[j]=1.0/(-2-V[j]+2.0i/s+alpha[j]);
        alpha[j-1]=-1.0*gamma[j];
    }

    beta[N-1]=0;
    xi[N]=xi[0]=0;
    for(l=0;l<pasos;l++)
    {
        for(j=0;j<=N;j++)
        {
            datos << j << ", " << abs(phi[j][l])*abs(phi[j][l]) << ", " << V[j] << endl;
        }
        datos << endl; 
        for(j=N-1;j>0;j--)
        {
            beta[j-1]=gamma[j]*(4i*phi[j][l]/s-beta[j]);
        }
        for(j=0;j<N;j++)
        {
            xi[j+1]=alpha[j]*xi[j]+beta[j];
        }
        for(j=1;j<N;j++)
        {
            phi[j][l+1]=xi[j]-phi[j][l];
        }
    }

    datos.close();
    return 0;   
}