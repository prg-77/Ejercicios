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
void actualizar_aceleracion (double r[][2], double a[][2], int n_part, double L);
void actualizar_posicion (double r[][2], double v[][2], double a[][2], double w[][2], double h, int n_part, double L);
void actualizar_velocidad (double v[][2], double w[][2], double a[][2], double h, int n_part);
double signo (double num);

int main()
{
    double r[20][2], v[20][2], a[20][2], w[20][2], h, t;
    double L;
    int i,j,k,n, n_part;
    ofstream datos;

    datos.open("datos_posicion.dat");

    t=0, h=0.002;
    n_part=20;
    L=10;

    inicializar(r,v,w,n_part,L);
    //escribir_en_fichero(r,v,n_part);

    actualizar_aceleracion(r,a,n_part,L);

    for(i=0;i<8000;i++)
  {
    if(n%200==0)
    {
        for(j=0;j<20 ; j++)
        {
        datos << r[j][0] << " , " << r[j][1] << endl;
        }
        datos << endl;
    }
    
    actualizar_posicion(r,v,a,w,h,n_part,L);
    actualizar_aceleracion(r,a,n_part,L); 
    actualizar_velocidad(v,w,a,h,n_part);
    
    t=t+h;
  }


    datos.close();
}

void inicializar (double r[][2], double v[][2], double w[][2], int n_part, double L)
{
    int i, j;
    double theta;

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

void actualizar_aceleracion (double r[][2], double a[][2], int n_part, double L)
{
    int i,j;
    double radio, ac, theta;
    double c, x, y;

    c=1e-2;
    
    for(i=0;i<n_part;i++)
    {
        a[i][0]=a[i][1]=0;

        for(j=0;j<n_part;j++)
        {
           ac=0;
           if(i!=j)
           {
                if(abs(r[i][0]-r[j][0])<L/2)
                  {
                    x=(r[j][0]-r[i][0]);
                    if(abs(r[i][1]-r[j][1])<L/2)
                    {
                        y=(r[j][1]-r[i][1]);
                    }
                    else
                    {
                        y=-signo(r[j][1]-r[i][1])*(L-abs(r[j][1]-r[i][1]));                        
                    }
                  }

                  else
                  {
                    x=-signo(r[j][0]-r[i][0])*(L-abs(r[i][0]-r[j][0]));
                    if(abs(r[i][1]-r[j][1])<L/2)
                    {
                        y=(r[j][1]-r[i][1]);                        
                    }
                    else
                    {
                        y=-signo(r[j][1]-r[i][1])*(L-abs(r[i][1]-r[j][1]));
                    }
                  }
                radio=sqrt(x*x+y*y);
                ac=-4.0*((pow(c,6)*6/pow(radio,7))-(pow(c,12)*12/pow(radio,13)));
                
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