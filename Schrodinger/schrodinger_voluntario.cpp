#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <complex>

#define N 500
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define n_ciclos 100
#define lambda 0.3
#define pasos 1500
#define experimentos 1000 

using namespace std;

double calcular_norma(complex<double> phi[N+1]);
void generar(double& k, double& s, double V[N+1], complex<double> phi[N+1], complex<double> gamma[N], complex<double> alpha[N]);
void calcular_beta(double s, complex<double> gamma[N], complex<double> beta[N], complex<double> phi[N+1], complex<double> alpha[N], complex<double> xi[N+1]);
double calcular_PD(complex<double> phi[N+1]);
int calcular_maximo(double PD[pasos]);
double calcular_x(complex<double> phi[N+1]);
double calcular_p(complex<double> phi[N+1]);
double calcular_V(double V[N+1], complex<double> phi[N+1]);
double calcular_error_x(complex<double> phi[N+1]);
double calcular_error_p(complex<double> phi[N+1]);
double calcular_error_V(double V[N+1], complex<double> phi[N+1]);
double calcular_T(complex<double> phi[N+1]);
double calcular_error_T(complex<double> phi[N+1]);
double calcular_E(double V[N+1], complex<double> phi[N+1]);
double calcular_error_E(double V[N+1], complex<double> phi[N+1]);

int main()
{

    double k, s, V[N+1], PD[pasos], p, merr;
    complex<double> phi[N+1], gamma[N], alpha[N], beta[N], xi[N+1];
    int j, l, m, mT, mi[experimentos];
    //ofstream datos;
    ofstream datos_observables, datos_PD;

    //datos.open("schrodinger_data.dat");
    datos_observables.open("graficas_observables.dat");
    datos_PD.open("PD.dat");
  
    mT=0;

    srand(time(NULL));

    for(m=0;m<experimentos;m++)
    {
        generar(k,s,V,phi,gamma,alpha);

        beta[N-1]=xi[N]=xi[0]=0;

        for(l=0;l<pasos;l++)
        {
            
            PD[l]=calcular_PD(phi);
            datos_PD << l << "  " << PD[l] << endl;

            //Crear datos para grafica de observables
            /*
            if(l%10==0)
            {
               datos_observables << l << "  " << calcular_x(phi) << "  " << calcular_error_x(phi) << "  " << calcular_p(phi) << "  " << calcular_error_p(phi) << "  " << calcular_V(V,phi) << "  " << calcular_error_V(V,phi) << "  " << calcular_T(phi) << "  " << calcular_error_T(phi) << "  " << calcular_E(V,phi) << "  " << calcular_error_E(V,phi) << endl;
            }
            */

            /*
            if(l%3==0)
            {
                for(j=0;j<=N;j++)
                {
                    datos << j << ", " << abs(phi[j])*abs(phi[j]) << ", " << V[j] << ", " << real(phi[j])*real(phi[j]) << ", " << imag(phi[j])*imag(phi[j]) << ", " <<  calcular_norma(phi) <<  endl;
                }
                datos << endl;
            }
            */
            
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
    merr=sqrt(merr/experimentos)/sqrt(experimentos);  // Intervalo de confianza del 68%

    cout << "Probabilidad de transmision: " << 1.0*mT/experimentos << " +/- " << merr << endl;

    //datos.close();
    datos_observables.close();
    datos_PD.close();
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
    double norma;
    int j;
    
    k=2.0*PI*n_ciclos/N;
    s=1.0/(4*k*k);

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
        phi[j]=pow(exp(1.0i),k*j)*exp(-8.0*(4*j-N)*(4*j-N)/(N*N));
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
    for(j=4*N/5;j<=N;j++)
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
            if(PD[j+2]>1.0001*PD[j+1] || PD[j+3]>1.0001*PD[j+2] || PD[j+4]>1.0001*PD[j+3] || PD[j+5]>1.0001*PD[j+4] || PD[j+6]>1.0001*PD[j+5])
            {
                j++;
            }
            else
            {
                mayor=false;
            }
            
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

// Cálculo del observable p (módulo)

double calcular_p(complex<double> phi[N+1])
{
    int j;
    complex<double> p;
    p=0;
    for(j=0;j<N;j++)
    {
        p=p+conj(phi[j])*(phi[j+1]-phi[j]);
    }
    return abs(p);
}

// Cálculo del observable T

double calcular_T(complex<double> phi[N+1])
{
    int j;
    complex<double> T_obs;
    T_obs=0;
    for(j=0;j<N;j++)
    {
        T_obs=T_obs-conj(phi[j])*(phi[j+2]-2.0*phi[j+1]+phi[j]);
    };
    return abs(T_obs);
}

// Cálculo del observable V

double calcular_V(double V[N+1], complex<double> phi[N+1])
{
    int j;
    double V_obs;
    V_obs=0;
    for(j=0;j<=N;j++)
    {
        V_obs=V_obs+V[j]*abs(phi[j])*abs(phi[j]);
    }
    return V_obs;
}

// Cálculo del observable E=T+V

double calcular_E(double V[N+1], complex<double> phi[N+1])
{
    int j;
    return calcular_T(phi)+calcular_V(V,phi);
}


//Error en el observable x (para suma de Riemann con la regla del punto medio)

double calcular_error_x(complex<double> phi[N+1])
{
    int j, maximo;
    double der2[N-2];
    for(j=0;j<(N-1);j++)
    {
        der2[j]=(j+2)*abs(phi[j+2])*abs(phi[j+2])-2*(j+1)*abs(phi[j+1])*abs(phi[j+1])+j*abs(phi[j])*abs(phi[j]);
    }
    maximo=0;
    for(j=1;j<(N-1);j++)
    {
        if(der2[j]>der2[maximo])
        {
            maximo=j;
        }
    }
    return der2[maximo]*N/24;
}

//Error en el observable p (para suma de Riemann con la regla del punto medio)

double calcular_error_p(complex<double> phi[N+1])
{
    int j, maximo;
    complex<double> der2[N-3];
    for(j=0;j<(N-2);j++)
    {
        der2[j]=conj(phi[j+2])*(phi[j+3]-phi[j+2])-2.0*conj(phi[j+1])*(phi[j+2]-phi[j+1])+conj(phi[j])*(phi[j+1]-phi[j]);
    }
    maximo=0;
    for(j=1;j<(N-2);j++)
    {
        if(abs(der2[j])>abs(der2[maximo]))
        {
            maximo=j;
        }
    }
    return abs(der2[maximo])*N/24;
}

//Error en el observable T (para suma de Riemann con la regla del punto medio)

double calcular_error_T(complex<double> phi[N+1])
{
    int j, maximo;
    complex<double> der2[N-3];
    for(j=0;j<(N-3);j++)
    {
        der2[j]=abs(conj(phi[j+2])*(phi[j+4]-2.0*phi[j+3]+phi[j+2]))-abs(2.0*conj(phi[j+1])*(phi[j+3]-2.0*phi[j+2]+phi[j+1]))+abs(conj(phi[j])*(phi[j+2]-2.0*phi[j+1]+phi[j]));
    }
    maximo=0;
    for(j=1;j<(N-3);j++)
    {
        if(abs(der2[j])>abs(der2[maximo]))
        {
            maximo=j;
        }
    }
    return abs(der2[maximo])*N/24;
}

//Error en el observable V (para suma de Riemann con la regla del punto medio)

double calcular_error_V(double V[N+1], complex<double> phi[N+1])
{
    int j, maximo;
    double der2[N-2];
    for(j=0;j<(N-1);j++)
    {
        der2[j]=V[j+2]*abs(phi[j+2])*abs(phi[j+2])-2*V[j+1]*abs(phi[j+1])*abs(phi[j+1])+V[j]*abs(phi[j])*abs(phi[j]);
    }
    maximo=0;
    for(j=1;j<(N-1);j++)
    {
        if(der2[j]>der2[maximo])
        {
            maximo=j;
        }
    }
    return der2[maximo]*N/24;
}

//Error en el observable E (para suma de Riemann con la regla del punto medio)

double calcular_error_E(double V[N+1], complex<double> phi[N+1])
{
    return calcular_error_T(phi)+calcular_error_V(V,phi);
}