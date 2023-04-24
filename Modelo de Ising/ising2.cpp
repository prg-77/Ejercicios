#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include "gsl_x64-windows/include/gsl/gsl_rng.h"

#define N 100
#define pasos 1000

using namespace std;

void inicializar(int s[][N+2]);
double calcular_deltaE (int s[][N+2], int n, int m);
double calcular_p (double deltaE, double T);
void aceptacion_rechazo (double p, int s[][N+2], int n, int m);
void generar_n_m (int& n, int& m);

gsl_rng*tau;

int main()
{
    extern gsl_rng*tau;

    double T,deltaE;
    double p;
    int i,j,k,n,m,s[N+2][N+2];
    ofstream datos;

    int semilla=18237247;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);

    datos.open("datos.dat");

    //T=0.001;
    T=2.3;


    inicializar(s);

    for(k=0;k<pasos;k++)
    {
          for(i=1;i<N+1;i++)
        {
            for(j=1;j<N;j++)
            {
                datos << s[i][j] << ", ";
            }
            datos << s[i][N] << endl;
        }

        datos << endl;
        
        for(i=0;i<N*N;i++)
        {    
            generar_n_m(n,m);

            deltaE=calcular_deltaE(s,n,m);

            p=calcular_p(deltaE,T);

            aceptacion_rechazo(p,s,n,m);
        }
    }
    datos.close();
}

void inicializar(int s[][N+2])
{
    int i,j;
    double t;
    
    srand(time(NULL));

    for(i=1;i<N+1;i++)
    {
        for(j=1;j<N+1;j++)
        {
            t=1.0*rand()/RAND_MAX;
            if(t<0.5)
                s[i][j]=1;
            else
                s[i][j]=-1;
            //s[i][j]=1;
        }
    }
}

double calcular_deltaE (int s[][N+2], int n, int m)
{
    if(n==1)
        s[n-1][m]=s[N][m];
    if(m==1)
        s[n][m-1]=s[n][N];
    if(n==N)
        s[n+1][m]=s[1][m];
    if(m==N)
        s[n][m+1]=s[n][1];

    return 2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
}

double calcular_p (double deltaE, double T)
{
    double p;

    if(1<exp(-deltaE/T))
    {
         p=1;
     }
    else
    {
        p=exp(-deltaE/T);
    }

    return p;
}

void aceptacion_rechazo (double p, int s[][N+2], int n, int m)
{
    double xi;

    xi=1.0*rand()/RAND_MAX;

    if(xi<p)
    {
        s[n][m]=-s[n][m];
    }

    return;
}

void generar_n_m (int& n, int& m)
{
    n=rand()%N+1;
    m=rand()%N+1;
}