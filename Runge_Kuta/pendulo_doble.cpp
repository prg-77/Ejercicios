#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define g 9.8
#define energia 10.0
#define pasos 2000
#define h 0.02

using namespace std;

double f_phi(double phi, double psi, double p_phi, double p_psi);
double f_p_phi(double phi, double psi, double p_phi, double p_psi);
double f_psi(double phi, double psi, double p_phi, double p_psi);
double f_p_psi(double phi, double psi, double p_phi, double p_psi);

int main()
{
    double y[3], t, k1[3], k2[3], k3[3], k4[3];
    int i, j;
    ofstream animacion, grafica;

    animacion.open("pendulo_doble.dat");
    grafica.open("grafica_pendulo.dat");

    //y[0]=phi, y[1]=p_phi, y[2]=psi, y[3]=p_psi

    // Inicializar variables

    t=0;
    y[0]=0.0; // Valor inicial de phi (en radianes)
    y[2]=0.0; // Valor inicial de psi (en radianes)
    y[1]=2*sqrt(energia-2*g*(1-cos(y[0]))-g*(1-cos(y[2]))); // Valor inicial de p_phi (en radianes)
    y[3]=sqrt(energia-2*g*(1-cos(y[0]))-g*(1-cos(y[2])))*cos(y[0]-y[2]);  // Valor inicial de  p_psi (en radianes)

    for(i=0;i<pasos;i++)
    {
       
       animacion << sin(y[0]) << ", " << -cos(y[0]) << endl;
       animacion << sin(y[0])+sin(y[2]) << ", " << -(cos(y[0])+cos(y[2])) << endl;
       animacion << endl; 
       grafica << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;

       k1[0]=h*f_phi(y[0],y[2],y[1],y[3]);
       k1[1]=h*f_p_phi(y[0],y[2],y[1],y[3]);
       k1[2]=h*f_psi(y[0],y[2],y[1],y[3]);
       k1[3]=h*f_p_psi(y[0],y[2],y[1],y[3]); 

       k2[0]=h*f_phi(y[0]+k1[0]/2.0,y[2]+k1[2]/2.0,y[1]+k1[1]/2.0,y[3]+k1[3]/2.0);
       k2[1]=h*f_p_phi(y[0]+k1[0]/2.0,y[2]+k1[2]/2.0,y[1]+k1[1]/2.0,y[3]+k1[3]/2.0);
       k2[2]=h*f_psi(y[0]+k1[0]/2.0,y[2]+k1[2]/2.0,y[1]+k1[1]/2.0,y[3]+k1[3]/2.0);
       k2[3]=h*f_p_psi(y[0]+k1[0]/2.0,y[2]+k1[2]/2.0,y[1]+k1[1]/2.0,y[3]+k1[3]/2.0);

       k3[0]=h*f_phi(y[0]+k2[0]/2.0,y[2]+k2[2]/2.0,y[1]+k2[1]/2.0,y[3]+k2[3]/2.0);
       k3[1]=h*f_p_phi(y[0]+k2[0]/2.0,y[2]+k2[2]/2.0,y[1]+k2[1]/2.0,y[3]+k2[3]/2.0);
       k3[2]=h*f_psi(y[0]+k2[0]/2.0,y[2]+k2[2]/2.0,y[1]+k2[1]/2.0,y[3]+k2[3]/2.0);
       k3[3]=h*f_p_psi(y[0]+k2[0]/2.0,y[2]+k2[2]/2.0,y[1]+k2[1]/2.0,y[3]+k2[3]/2.0);

       k4[0]=h*f_phi(y[0]+k3[0],y[2]+k3[2],y[1]+k3[1],y[3]+k3[3]);
       k4[1]=h*f_p_phi(y[0]+k3[0],y[2]+k3[2],y[1]+k3[1],y[3]+k3[3]);
       k4[2]=h*f_psi(y[0]+k3[0],y[2]+k3[2],y[1]+k3[1],y[3]+k3[3]);
       k4[3]=h*f_p_psi(y[0]+k3[0],y[2]+k3[2],y[1]+k3[1],y[3]+k3[3]);

       for(j=0;j<4;j++)
       {
           y[j]=y[j]+(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/6.0;
       }

        t=t+h;
    }

    animacion.close();
    grafica.close();

    return 0;
}


double f_phi(double phi, double psi, double p_phi, double p_psi)
{
    return (p_phi-p_psi*cos(phi-psi))/(1.0+sin(phi-psi)*sin(phi-psi));
}

double f_p_phi(double phi, double psi, double p_phi, double p_psi)
{
    return -2*g*sin(phi)-(p_phi*p_psi*sin(phi-psi))/(1.0+sin(phi-psi)*sin(phi-psi))+sin(2*(phi-psi))*(p_phi*p_phi+2*p_psi*p_psi-2*p_phi*p_psi*cos(phi-psi))/(2.0*(1.0+sin(phi-psi)*sin(phi-psi))*(1.0+sin(phi-psi)*sin(phi-psi)));
}

double f_psi(double phi, double psi, double p_phi, double p_psi)
{
    return (2*p_psi-p_phi*cos(phi-psi))/(1.0+sin(phi-psi)*sin(phi-psi));
}

double f_p_psi(double phi, double psi, double p_phi, double p_psi)
{
    return -g*sin(psi)+(p_phi*p_psi*sin(phi-psi))/(1.0+sin(phi-psi)*sin(phi-psi))-sin(2*(phi-psi))*(p_phi*p_phi+2*p_psi*p_psi-2*p_phi*p_psi*cos(phi-psi))/(2.0*(1.0+sin(phi-psi)*sin(phi-psi))*(1.0+sin(phi-psi)*sin(phi-psi)));
}