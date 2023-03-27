#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

#define c 1.496e+11
#define Ms 1.989e+30
#define G 6.67e-11 

using namespace std;

void inicializar (double r[][2], double v[][2], double w[][2], double m[]);
void leerficheroyreescalar (double r[][2], double v[][2], double m[]);
void escribirreescaladofichero(double r[][2], double v[][2], double m[]);
void actualizar_aceleracion (double r[][2], double a[][2], double m[]);
void actualizar_posicion (double r[][2], double v[][2], double a[][2], double w[][2], double h);
void actualizar_velocidad (double v[][2], double w[][2], double a[][2], double h);

int main()
{
  double r[9][2], v[9][2], a[9][2], w[9][2], m[9],  h, t;
  int i,j,k,n;  
  ofstream datos;

  h=0.1, t=0;

  inicializar(r,v,w,m);
  leerficheroyreescalar(r,v,m);
  escribirreescaladofichero(r,v,m);
  
  actualizar_aceleracion(r,a,m);

  datos.open("datos_posiciones_experimental.dat");
  
  n=0; 

  for(i=0;i<10000;i++)
  {
    if(n%10==0)
    {
        for(j=0;j<9 ; j++)
        {
        datos << r[j][0]-r[2][0] << " , " << r[j][1]-r[2][1] << endl;
        }
        datos << endl;
    } 

    actualizar_posicion(r,v,a,w,h);
    actualizar_aceleracion(r,a,m); 
    actualizar_velocidad(v,w,a,h);
    
    t=t+h;
    n=n+1;
  }

  datos.close();

  return 0;  
}


void inicializar (double r[][2], double v[][2], double w[][2], double m[])
{
    int i, j;
    for(i=0;i<9;i++)
    {
        m[i]=0;
        for(j=0;j<2;j++)
            r[i][j]=v[i][j]=w[i][j]=0;
    }
    return;
}

void leerficheroyreescalar (double r[][2], double v[][2], double m[])
{
    int i;
    ifstream inicial;
    inicial.open("condicionesexperimentales.txt");
    for(i=0; i<9; i++)
  {
        inicial >> r[i][0];
        r[i][0]=r[i][0]/c;
        inicial >> v[i][1];
        v[i][1]=v[i][1]*(1/c)/pow(G*Ms/pow(c,3),0.5);
        inicial >> m[i];
        m[i]=m[i]/Ms;
  }
    //r[8][1]=r[8][0];
    //v[8][0]=v[8][1]+0.1e4*(1/c)/pow(G*Ms/pow(c,3),0.5);
    inicial.close();
    return;
}

void escribirreescaladofichero(double r[][2], double v[][2], double m[])
{
    int i;
    ofstream reescalado;

    reescalado.open("reescaladoexperimental.txt");
    
    for(i=0;i<9;i++)
    {
        reescalado << r[i][0] << "   ";
        reescalado << v[i][1] << "  ";
        reescalado << m[i] << "  " << endl;
    }

    reescalado.close();
    return;
}

void actualizar_aceleracion (double r[][2], double a[][2], double m[])
{
    int i,j,k;
    
    for(i=0;i<9;i++)
    {
        a[i][0]=a[i][1]=0;

        for(j=0;j<9;j++)
        {
           if(i!=j)
           {
                for(k=0;k<2;k++)
                {
                  a[i][k]=a[i][k]-m[j]*(r[i][k]-r[j][k])/pow(pow(r[i][0]-r[j][0],2)+pow(r[i][1]-r[j][1],2),1.5); 
                }
           }
        }
    }
    
    return;
}

void actualizar_posicion (double r[][2], double v[][2], double a[][2], double w[][2], double h)
{
    int i, j;

    for(i=0;i<9;i++)
    {
        for(j=0;j<2;j++)
        {
          w[i][j]=v[i][j]+h*a[i][j]/2; 
          r[i][j]=r[i][j]+h*w[i][j]; 
        }
    }
    return;
}

void actualizar_velocidad (double v[][2], double w[][2], double a[][2], double h)
{
    int i,j;

    for(i=0;i<9;i++)
    {
        for(j=0;j<2;j++)
        {
            v[i][j]=w[i][j]+h*a[i][j]/2;
        }
    }

    return;
}