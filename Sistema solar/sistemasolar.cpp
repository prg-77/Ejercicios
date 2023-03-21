#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

#define c 1.496e+11
#define Ms 1.989e+30
#define G 6.67e-11 

using namespace std;

void inicializar (double r[][2], double v[][2], double m[]);
void leerdatosiniciales (double r[][2], double v[][2], double m[]);
void escribirreescaladofichero(double r[][2], double v[][2], double m[]);

int main()
{
  double r[8][2], v[8][2], a[8][2], w[8][2], m[8],  h;

  inicializar(r,v,m);
  leerdatosiniciales(r,v,m);
  escribirreescaladofichero(r,v,m);

  return 0;
}


void leerdatosiniciales (double r[][2], double v[][2], double m[])
{
    int i;
    ifstream inicial;
    inicial.open("condicionesiniciales.txt");
    for(i=0; i<8; i++)
  {
        inicial >> r[i][0];
        r[i][0]=r[i][0]/c;
        inicial >> v[i][1];
        v[i][1]=v[i][1]*(1/c)/pow(G*Ms/pow(c,3),1/2);
        inicial >> m[i];
        m[i]=m[i]/Ms;
  }
    inicial.close();
    return;
}

void inicializar (double r[][2], double v[][2], double m[])
{
    int i, j;
    for(i=0;i<8;i++)
    {
        m[i]=0;
        for(j=0;j<2;j++)
            r[i][j]=v[i][j]=0;
    }
    return;
}

void escribirreescaladofichero(double r[][2], double v[][2], double m[])
{
    int i;
    ofstream reescalado;

    reescalado.open("reescalado.txt");
    
    for(i=0;i<8;i++)
    {
        reescalado << r[i][0] << "   ";
        reescalado << v[i][1] << "  ";
        reescalado << m[i] << "  " << endl;
    }

    reescalado.close();
    return;
}