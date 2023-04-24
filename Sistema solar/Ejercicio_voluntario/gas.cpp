#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265

using namespace std;

void inicializar (double r[][2], double v[][2], double w[][2], int n_part, double L);
void escribir_en_fichero(double r[][2], double v[][2], int n_part);
void actualizar_aceleracion (double r[][2], double a[][2], int n_part, double L, double sigma);
void actualizar_posicion (double r[][2], double v[][2], double a[][2], double w[][2], double h, int n_part, double L);
void actualizar_velocidad (double v[][2], double w[][2], double a[][2], double h, int n_part);
double signo (double num);
double calcular_radio (double r1[2], double r2[2], double L);
void energia (double& T, double& V, double n_part, double r[][2], double v[][2], double sigma, double L);
double temperatura (double T, double n_part);


int main()
{
    double r[20][2], v[20][2], a[20][2], w[20][2], h, t, sigma;
    double L, T, V, T_med;
    int i,j,k,n, n_part;
    ofstream datos, datos_energia;

    datos.open("datos_posicion.dat");
    datos_energia.open("datos_energia.dat");

    //sigma=1e-2;
    sigma=1e-2; //Debería ser 1 
    t=0, h=0.001;
    n_part=20;
    L=10;
    n=0;
    T_med=0;

    inicializar(r,v,w,n_part,L);
    //escribir_en_fichero(r,v,n_part);

    actualizar_aceleracion(r,a,n_part,L,sigma);

    for(i=0;i<5000;i++)
  {
    if(n%5==0)
    {
        for(j=0;j<20 ; j++)
        {
        datos << r[j][0] << " , " << r[j][1] << endl;
        }
        datos << endl;
        energia(T,V,n_part,r,v,sigma,L);
        datos_energia << T << "  " << V << "  " << T+V << endl;
        T_med=T_med+temperatura(T,n_part);
    }
    
    actualizar_posicion(r,v,a,w,h,n_part,L);
    actualizar_aceleracion(r,a,n_part,L,sigma); 
    actualizar_velocidad(v,w,a,h,n_part);
    
    t=t+h;
    n=n+1;
  }

    T_med=T_med/(n/5); // Cuidado al cambiar n% !!
    cout << "Temperatura media: " << T_med << endl;

    datos.close();
    datos_energia.close();
}

void inicializar (double r[][2], double v[][2], double w[][2], int n_part, double L)
{
    int i, j,k,l,n;
    double theta;
    bool separado;

    srand(time(NULL));

    for(i=0;i<n_part;i++)
    {
        theta=2*PI*rand()/RAND_MAX;
        v[i][0]=cos(theta);
        v[i][1]=sin(theta);
        for(j=0;j<2;j++)
        {
            r[i][j]=L*rand()/RAND_MAX;
            w[i][j]=0;
        }
        separado=false;
        while(separado==false)
        {
            separado=true;
            for(k=0;k<i;k++)
            {
                if(n>1000)
                {
                    i=0;
                    n=0;
                }
                if(calcular_radio(r[i],r[k],L)<=2)
                {
                    separado=false;
                    r[i][0]=L*rand()/RAND_MAX;
                    r[i][1]=L*rand()/RAND_MAX;
                    k=i;
                    n=n+1;
                }
            }
        }
            
    }
    return;
}

void escribir_en_fichero(double r[][2], double v[][2], int n_part)
{
    int i;
    ofstream fich;

    fich.open("datos_iniciales.txt");
    
    for(i=0;i<n_part;i++)
    {
        fich << r[i][0] << "  " << r[i][1] << "  ";
        fich << v[i][0] << "  " << v[i][1] << "  ";
        fich << endl;
    }

    fich.close();
    return;
}

void actualizar_aceleracion (double r[][2], double a[][2], int n_part, double L, double sigma)
{
    int i,j;
    double radio, ac, theta;
    double x, y;
    
    for(i=0;i<n_part;i++)
    {
        a[i][0]=a[i][1]=0;

        for(j=0;j<n_part;j++)
        {
           ac=0;
           if(i!=j)
           {
                radio=calcular_radio(r[i],r[j],L);
                ac=-4.0*((pow(sigma,6)*6/pow(radio,7))-(pow(sigma,12)*12/pow(radio,13)));
                
                if(x!=0)
                {
                    if(x>0) 
                    {
                        theta=atan(y/x);
                    }
                    else
                    {
                        theta=atan(y/x)+PI;
                    }
                }
                else
                {
                    if(y>0)
                        theta=PI/2;
                    else if (y<0)
                        theta=-PI/2;
                    else 
                        ac=0;
                }
                
                a[i][0]=a[i][0]+ac*cos(theta);
                a[i][1]=a[i][1]+ac*sin(theta);
                }
           }
        }

    
    return;
}

void actualizar_posicion (double r[][2], double v[][2], double a[][2], double w[][2], double h, int n_part, double L)
{
    int i, j;
    double b;

    for(i=0;i<n_part;i++)
    {
        for(j=0;j<2;j++)
        {
          w[i][j]=v[i][j]+h*a[i][j]/2; 
          r[i][j]=r[i][j]+h*w[i][j];

          if(r[i][j]>L)
          {
           b=r[i][j]/L;
           r[i][j]=r[i][j]-L*b;
          }      
           
          else if(r[i][j]<0)
          {
                b=-r[i][j]/L;
                r[i][j]=r[i][j]+(b+1)*L;
          }
    
        }
    }

    return;
}

void actualizar_velocidad (double v[][2], double w[][2], double a[][2], double h, int n_part)
{
    int i,j;

    for(i=0;i<n_part;i++)
    {
        for(j=0;j<2;j++)
        {
            v[i][j]=w[i][j]+h*a[i][j]/2;
        }
    }

    return;
}

double signo (double num)
{
    if(num>0)
        return 1.0;
    else if(num<0)
        return -1.0;
    else
        return 0.0;
}

double calcular_radio (double r1[2], double r2[2], double L)
{
    double x, y;
    if(abs(r1[0]-r2[0])<L/2)
    {
            x=(r2[0]-r1[0]);
            if(abs(r2[1]-r1[1])<L/2)
            {
                y=(r2[1]-r1[1]);
            }
            else
            {
                y=-signo(r2[1]-r1[1])*(L-abs(r2[1]-r1[1]));                        
            }
            }

            else
            {
            x=-signo(r2[0]-r1[0])*(L-abs(r1[0]-r2[0]));
            if(abs(r1[1]-r2[1])<L/2)
            {
                y=(r2[1]-r1[1]);                        
            }
            else
            {
                y=-signo(r2[1]-r1[1])*(L-abs(r1[1]-r2[1]));
            }
    }
        return sqrt(x*x+y*y);                
}

void energia (double& T, double& V, double n_part, double r[][2], double v[][2], double sigma, double L)
{
    int i,j;
    double radio;

    T=V=0;
    for(i=0;i<n_part;i++)
    {
        T=T+(v[i][0]*v[i][0]+v[i][1]*v[i][1])/2;
        for(j=0;j<n_part;j++)
        {
            if(i!=j)
            {
                radio=calcular_radio(r[i],r[j],L);
                V=V+4*((pow(sigma,12)/pow(radio,12))-(pow(sigma,6)/pow(radio,6)));
            }
        }
    }
    return;
}

double temperatura (double T, double n_part)
{
    return T/n_part;
}