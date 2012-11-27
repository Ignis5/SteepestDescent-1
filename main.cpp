#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <iomanip>
#define N 2
using namespace std;

// funkcija six hump
double sixhump(double * x)
{
    double f = (4- 2.1 * x[0] * x[0] + (pow(x[0],4))/3) * x[0] * x[0] + x[0] * x[1] + ((-4+4*x[1]*x[1])*(pow(x[1],2)));
    return f;
}
double sixhump(double x, double y)
{
    double f = (4- 2.1 * x * x + (pow(x,4))/3) * x * x + x * y + ((-4+4*y*y)*(pow(y,2)));
    return f;
}

struct taskai{
    double x;
    double y;
    double f;
};

// Vektoriaus begalines (max) normos funkcijos deklaracija
double Vector_Max_Norm(double v[], int n);

// Greiciausio nusileidimo (angl. Steepest Descent) metodo deklaracija
int  Steepest_Descent(double (*f)(double *), void (*df)(double *, double *),
     int (*stopping_rule)(double*, double, double*, double, double*, int, int),
                          double a[], double *fa, double *dfa, double cutoff,
						double cutoff_scale_factor, double tolerance, int n);

// Generuoja atsitiktini realu skaiciu tarp dLow and dHigh
double GetRandomNumber(double dLow, double dHigh){
    return static_cast<double>(rand())/RAND_MAX*(dHigh-dLow) + dLow;
}

// Apskaiciuoja Six-hump Camel Back gradiento reiksme taske x
void SixHumpCamelBackGradient(double *x, double *fGrad){
    fGrad[0] = 8*x[0]-8.4*x[0]*x[0]*x[0]+2*x[0]*x[0]*x[0]*x[0]*x[0]+x[1];
    fGrad[1] = x[0]-8*x[1]+16*x[1]*x[1]*x[1];
}

// Algoritmo sustojimo salyga kontroliuojanti funkcija
int StoppingRule(double* a, double fa, double* x, double fx, double* dfa, int
iteration, int n){
	double fEps = abs(fx - fa); // Funkcijos reiksmiu skirtumas
	double xa[n];
	for(int i = 0; i < n; ++i) xa[i] = x[i]-a[i];
	double xEps = Vector_Max_Norm(xa, 2); // Argumento skirtumo norma
	double dfaEps = Vector_Max_Norm(dfa, 2); // Gradiento norma
	if(iteration > 3)
		return -6;
	else
		return 0;
}


int main(){

    srand(time(0)); // Kad kiekviena karta randomas suveiktu skirtingai
    cout << "Monte Carlo realizacijos pradzia" << endl;
    cout << "(RANDOM SEARCH METHOD)" << endl;
    // Kintamuju apsirasymas
    double x1 = -1.9, y1 = -1.1, x2 = 1.9, y2 = 1.1; // apsirasome intervalo rezius
    double EPS = 0.01, sprendinys = 0;
    int i, j, k;
    vector <taskai> r(0);
    struct taskai laik;
    double max;
    double tikrasis_sprendinys = -1.031628453;
    int indx;
    int kiek;
    cout << "Kiek surasti artimiausiu sprendiniui reiksmiu? ";
    cin >> kiek;
    double min[kiek][3];
    double a[N];

    j = 0;
    while (abs(-1.031628453 - sprendinys) > EPS){
        j = j + 1;

        laik.x = rand()*(x2 - x1)/(double)RAND_MAX+x1;
        laik.y = rand()*(y2 - y1)/(double)RAND_MAX+y1;
        a[0] = laik.x;
        a[1] = laik.y;
        sprendinys = sixhump(a);
        laik.f = sprendinys;
        r.push_back(laik);
    }

    cout << "sprendinys sustabdes cikla: " << endl << fixed << setprecision(9) <<  sprendinys << endl;


    //surandame tolimiausia reiksme;
    max = r[1].f;
    for(int i=1; i<=j+1; i++)
	{
		if(max < r[i].f)
		{
			max = r[i].f;
		}
	}
	//cout << max << endl;
    for (i=1; i<=kiek; i++){
        min[i-1][0] = r[0].x;
        min[i-1][1] = r[0].y;
        min[i-1][2] = r[0].f;
        indx = 0;
        for (k=2; k< j+1; k++)
        {
            if (abs(-1.031628453 - min[i-1][2]) > abs(-1.031628453 - r[k].f)){
                min[i-1][0] = r[k].x;
                min[i-1][1] = r[k].y;
                min[i-1][2] = r[k].f;
                indx = k;
            }
        }
        r[indx].f = max;
    }
    cout << "Tikrasis sprendinys " << endl << tikrasis_sprendinys << endl;
    cout << "Artimiausios reiksmes: " << endl;
    for (i=0;i<kiek;i++)
        cout << "[" << min[i][0] << " " << min[i][1] << "] f = " <<  min[i][2] << endl;

    x1 = min[0][0];
    x2 = min[0][0];
    y1 = min[0][1];
    y2 = min[0][1];

    for (i=1;i<kiek;i++){
        if (x1 > min[i][0]) x1 = min[i][0];
        else if (x2 < min[i][0]) x2 = min[i][0];
        if (y1 > min[i][1]) y1 = min[i][1];
        else if (y2 < min[i][1]) y2 = min[i][1];
    }




    double region[] = {-1.9, 1.9, -1.1, 1.1};
	// N-matis Vektorius
	/*srand(time(0)); // Naudoja vis kita seed'a
	double a[N]; // N-matis Vektorius
	for(int i = 0; i < N; ++i){
        a[i] = GetRandomNumber(region[2*i], region[2*i+1]);
    }*/

for (i = 0; i < kiek; i++){

    a[0] = min[i][0];
    a[1] = min[i][1];

    cout << "Pateikiama " << a[0] << " " << a[1] << endl;

	double fa = min[i][2]; // Funkcijos reiksme pradiniame taske a
	double dfa[N];
	SixHumpCamelBackGradient(a, dfa); // Funkcijos gradiento reiksme taske a
	double cutoff = 1.0, cutoff_scale_factor = 1.0; // Pap. parametrai
	double tolerance = 0.01;
	int err = Steepest_Descent( sixhump, SixHumpCamelBackGradient, StoppingRule, a, &fa, dfa, cutoff, cutoff_scale_factor, tolerance, N);
	cout << endl;
	switch (err)
	{
		case 0:
			cout << "Success" << endl;
			break;
		case -1:
			cout << "In the line search three points are collinear." << endl;
			break;
		case -2:
			cout << "In the line search the extremum of the parabola through the three points is a maximum." << endl;
			break;
		case -3:
			cout << "Int the line search the initial points failed to satisfy the condition that x1 < x2 < x3 and fx1 > fx2 < fx3." << endl;
			break;
		case -4:
			cout << "Not enough HEAP memory." << endl;
			break;
		case -5:
			cout << "The gradient evaluated at the initial point vanishes." << endl;
		case -6:
			cout << "Exceed maximal number of iterations." << endl;
		break;
	}
	cout << "Greiciausio nusileidimo (angl. Steepest Descent) metodu" << endl;
	cout << "surastas sprendinys yra:" << endl;
	cout << "xMin = (" << a[0] << ", " << a[1] << ")" << endl;
	cout << "f(xMin) = " << fa << endl << endl;
	cout << "==============================================" << endl << endl;

}
return 0;
}
