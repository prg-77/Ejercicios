#include <iostream>  // contiene cin y cout
#include <vector>
#include <string> // para trabajar con cadenas de caracteres
#include <cmath> // contiene funciones matematicas
#include <fstream> //contiene funciones para tratamiento de ficheros

#define PI 3.1415

using namespace std;

int suma(int x, int y);

int main()
{
    int i, x, y;
    
    cout << tan(7*PI/6) << endl;
    cout << tan(-5*PI/6) << endl;
    /*
    cout << "Hello" << endl; 
    cout << "Introduce x: ";
    cin >> x;
    cout << "Introduce y: ";
    cin >> y;
    cout << "La suma es: " << suma(x,y) << endl;
    for(i=1; i<=50; i++)
    {
        cout << i << endl;
    }  
    */
    return 0;
}

int suma (int x, int y)
{
    return x+y;
}
