#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define G 6.67e-11
#define M_T 5.9736e24
#define M_L 7.349e22
#define d_TL 3.844e8
#define omega 2.6617e-6
#define R_T 6.37816e6
#define R_L 1.7374e6

using namespace std;

double f_r(double p_r);
double f_p_r(double r, double phi, double p_phi, double t, double delta, double mu);
double f_phi(double r, double p_phi);
double f_p_phi(double r, double phi, double t, double delta, double mu);
double calcular_r_prima(double r, double phi, double t);
double calcular_H (double r, double p_r, double phi, double p_phi, double t);

int main()
{

double delta, mu, y[3], masa, h, t, k1[3], k2[3], k3[3], k4[3];
int i, j, pasos;

//y[0]=r, y[1]=p_r, y[2]=phi, y[3]=p_phi

delta=G*M_T/pow(d_TL,3);
mu=M_L/M_T;
ofstream datos, energia;

//Inicializar variables

pasos=10000;
t=0;
h=60; 
masa=2.97e+6; 
y[0]=R_T/d_TL;
y[1]=11.092e+3/d_TL; //Cambiar 1 por momento inicial
y[2]=47.0*PI/180; //Cambiar por angulo inicial
y[3]=0.0/(masa*pow(d_TL,2)); //Cambiar 1 por momento inicial


datos.open("rungekuta.dat");
energia.open("energia.dat");

for(i=0;i<pasos;i++)
{
    if(i%20==0)
    {
        datos << 0 << ", " << 0 << endl;
        datos << y[0]*cos(y[2]) << ", " << y[0]*sin(y[2]) << endl;
        datos << 1.0*cos(omega*t) << ", " << 1.0*sin(omega*t) << endl;
        datos << endl;
        energia << calcular_H(y[0], y[1], y[2], y[3], t) << endl;
    }
    
    k1[0]=h*f_r(y[1]);
    k1[1]=h*f_p_r(y[0], y[2], y[3], t, delta, mu);
    k1[2]=h*f_phi(y[0], y[3]);
    k1[3]=h*f_p_phi(y[0], y[2], t, delta, mu);

    k2[0]=h*f_r(y[1]+k1[1]/2);
    k2[1]=h*f_p_r(y[0]+k1[0]/2, y[2]+k1[2]/2, y[3]+k1[3]/2, t+h/2, delta, mu);
    k2[2]=h*f_phi(y[0]+k1[0]/2, y[3]+k1[3]/2);
    k2[3]=h*f_p_phi(y[0]+k1[0]/2, y[2]+k1[2]/2, t+h/2, delta, mu);

    k3[0]=h*f_r(y[1]+k2[1]/2);
    k3[1]=h*f_p_r(y[0]+k2[0]/2, y[2]+k2[2]/2, y[3]+k2[3]/2, t+h/2, delta, mu);
    k3[2]=h*f_phi(y[0]+k2[0]/2, y[3]+k2[3]/2);
    k3[3]=h*f_p_phi(y[0]+k2[0]/2, y[2]+k2[2]/2, t+h/2, delta, mu);

    k4[0]=h*f_r(y[1]+k3[1]);
    k4[1]=h*f_p_r(y[0]+k3[0], y[2]+k3[2], y[3]+k3[3], t+h, delta, mu);
    k4[2]=h*f_phi(y[0]+k3[0], y[3]+k3[3]);
    k4[3]=h*f_p_phi(y[0]+k3[0], y[2]+k3[2], t+h, delta, mu);

    for(j=0;j<4;j++)
    {
        y[j]=y[j]+(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;
    }

    t=t+h;

}

datos.close();
energia.close();

    return 0;
}

double f_r(double p_r)
{
    return p_r;
}

double f_p_r(double r, double phi, double p_phi, double t, double delta, double mu)
{
    double r_prima;
    r_prima=calcular_r_prima(r, phi, t);
    return (p_phi*p_phi/pow(r,3))-delta*((1.0/pow(r,2))+(mu/pow(r_prima,3))*(r-cos(phi-omega*t)));
}

double f_phi(double r, double p_phi)
{
    return p_phi/(r*r);
}

double f_p_phi(double r, double phi, double t, double delta, double mu)
{
    double r_prima;
    r_prima=calcular_r_prima(r, phi, t);
    return -delta*r*(mu/pow(r_prima,3))*sin(phi-omega*t);
}

double calcular_r_prima(double r, double phi, double t)
{
    return sqrt(1+pow(r,2)-2*r*cos(phi-omega*t));
}

double calcular_H (double r, double p_r, double phi, double p_phi, double t)
{   
    double r_L;

    r_L=calcular_r_prima(r, phi, t);
    
    return (p_r*p_r/(2))+(p_phi*p_phi/(2*r*r))-(1/r)-(1/r_L)-omega*p_phi;
}