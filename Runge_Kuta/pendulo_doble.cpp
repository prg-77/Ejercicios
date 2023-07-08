#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define g 9.8
#define energia 23
#define pasos 2000
#define h 0.01

using namespace std;

double f_phi(double phi, double psi, double p_phi, double p_psi);
double f_p_phi(double phi, double psi, double p_phi, double p_psi);
double f_psi(double phi, double psi, double p_phi, double p_psi);
double f_p_psi(double phi, double psi, double p_phi, double p_psi);
void runge_kutta(double y[4]);
double calcular_phi_prima(double y[4]);
double calcular_psi_prima(double y[4]);

int main()
{
    double k1[4],k2[4],k3[4],k4[4];
    double y[4], t, d_0[2], lambda[2][pasos], y2[4], lya_phi, lya_psi, err_lya_phi, err_lya_psi;
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

    y2[0]=y[0]+0.0000001;
    y2[2]=y[2]+0.0000001;
    y2[1]=2*sqrt(energia-2*g*(1-cos(y2[0]))-g*(1-cos(y2[2])));
    y2[3]=sqrt(energia-2*g*(1-cos(y2[0]))-g*(1-cos(y2[2])))*cos(y2[0]-y2[2]);
    
    d_0[0]=sqrt(pow(y2[0]-y[0],2)+pow(calcular_phi_prima(y2)-calcular_phi_prima(y),2));
    d_0[1]=sqrt(pow(y2[2]-y[2],2)+pow(calcular_psi_prima(y2)-calcular_psi_prima(y),2));
    lya_phi=lya_psi=0.0;

    for(i=0;i<pasos;i++)
    {
       
    animacion << sin(y[0]) << ", " << -cos(y[0]) << endl;
    animacion << sin(y[0])+sin(y[2]) << ", " << -(cos(y[0])+cos(y[2])) << endl;
    animacion << endl; 
    grafica << y[0] << " " << calcular_phi_prima(y) << " " << y[2] << " " << calcular_psi_prima(y) << endl;


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

    k1[0]=h*f_phi(y2[0],y2[2],y2[1],y2[3]);
    k1[1]=h*f_p_phi(y2[0],y2[2],y2[1],y2[3]);
    k1[2]=h*f_psi(y2[0],y2[2],y2[1],y2[3]);
    k1[3]=h*f_p_psi(y2[0],y2[2],y2[1],y2[3]);

    k2[0]=h*f_phi(y2[0]+k1[0]/2.0,y2[2]+k1[2]/2.0,y2[1]+k1[1]/2.0,y2[3]+k1[3]/2.0);
    k2[1]=h*f_p_phi(y2[0]+k1[0]/2.0,y2[2]+k1[2]/2.0,y2[1]+k1[1]/2.0,y2[3]+k1[3]/2.0);
    k2[2]=h*f_psi(y2[0]+k1[0]/2.0,y2[2]+k1[2]/2.0,y2[1]+k1[1]/2.0,y2[3]+k1[3]/2.0);
    k2[3]=h*f_p_psi(y2[0]+k1[0]/2.0,y2[2]+k1[2]/2.0,y2[1]+k1[1]/2.0,y2[3]+k1[3]/2.0);

    k3[0]=h*f_phi(y2[0]+k2[0]/2.0,y2[2]+k2[2]/2.0,y2[1]+k2[1]/2.0,y2[3]+k2[3]/2.0);
    k3[1]=h*f_p_phi(y2[0]+k2[0]/2.0,y2[2]+k2[2]/2.0,y2[1]+k2[1]/2.0,y2[3]+k2[3]/2.0);
    k3[2]=h*f_psi(y2[0]+k2[0]/2.0,y2[2]+k2[2]/2.0,y2[1]+k2[1]/2.0,y2[3]+k2[3]/2.0);
    k3[3]=h*f_p_psi(y2[0]+k2[0]/2.0,y2[2]+k2[2]/2.0,y2[1]+k2[1]/2.0,y2[3]+k2[3]/2.0);

    k4[0]=h*f_phi(y2[0]+k3[0],y2[2]+k3[2],y2[1]+k3[1],y2[3]+k3[3]);
    k4[1]=h*f_p_phi(y2[0]+k3[0],y2[2]+k3[2],y2[1]+k3[1],y2[3]+k3[3]);
    k4[2]=h*f_psi(y2[0]+k3[0],y2[2]+k3[2],y2[1]+k3[1],y2[3]+k3[3]);
    k4[3]=h*f_p_psi(y2[0]+k3[0],y2[2]+k3[2],y2[1]+k3[1],y2[3]+k3[3]);

    for(j=0;j<4;j++)
    {
        y2[j]=y2[j]+(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/6.0;
    }
       
    //Coeficientes de Lyapunov

    lambda[0][i]=log(sqrt(pow(y2[0]-y[0],2)+pow(calcular_phi_prima(y2)-calcular_phi_prima(y),2))/d_0[0]);
    lambda[1][i]=log(sqrt(pow(y2[2]-y[2],2)+pow(calcular_psi_prima(y2)-calcular_psi_prima(y),2))/d_0[1]);

    lya_phi=lya_phi+lambda[0][i];
    lya_psi=lya_psi+lambda[1][i];   

        t=t+h;
    }

    lya_phi=lya_phi/t;
    lya_psi=lya_psi/t;

    err_lya_phi=err_lya_psi=0.0;
    for(i=0;i<pasos;i++)
    {
        err_lya_phi=err_lya_phi+pow(lambda[0][i]-lya_phi,2);
        err_lya_psi=err_lya_psi+pow(lambda[1][i]-lya_psi,2);
    }

    err_lya_phi=sqrt(err_lya_phi/pasos)/sqrt(pasos);
    err_lya_psi=sqrt(err_lya_psi/pasos)/sqrt(pasos);

    animacion.close();
    grafica.close();
    cout << "Coeficiente de Lyapunov phi: " << lya_phi << " +/- " << err_lya_phi << endl;
    cout << "Coeficiente de Lyapunov psi: " << lya_psi << " +/- " << err_lya_psi << endl;

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

void runge_kutta(double y[4])
{
    double k1[3],k2[3],k3[3],k4[3];
    int j;

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

    return; 
}

double calcular_phi_prima(double y[4])
{
    return (y[1]-y[3]*cos(y[0]-y[2]))/(1.0+sin(y[0]-y[2])*sin(y[0]-y[2]));
}

double calcular_psi_prima(double y[4])
{
    return (2.0*y[3]-y[1]*cos(y[0]-y[2]))/(1.0+sin(y[0]-y[2])*sin(y[0]-y[2]));
}
