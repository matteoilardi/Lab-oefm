#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>

#define N_SETS 3
#define C 299880000
#define N_P 10000000
#define _USE_MATH_DEFINES

using namespace std;


struct measure{
    double nu;
    double omega;
    double delta;
};

struct set{
    measure* v;
    int dim;
};

struct result{
    double m;
    double stdev;

    result();
    result(double* v, int n);
    result(result* vr, int n);
};

void ignoreText(ifstream* is);
void readMeasure(ifstream* is, measure* m);
void printMeasure(ofstream* os, const measure m);
void printSet(ofstream *os, const set s);
void reldif(const set s, double* omega, double* delta, int n); //sottrae l'n-esimo valore al primo
void reldif2(const set s, double* omega, double* delta, int n); //sottrae i valori a coppie di due, usata per il terzo set

void setResult(result& r, double* v, int n);

double dc_da(double c, double F, double a, double d);
double dc_dd(double c, double F, double a, double d);
double dc_dF(double c, double F, double a, double d);
double domega_ddelta(double c, double F, double a, double d);

double rand_double(double min, double max);
double stGaussian(double z);
double stGauusianIntegral(int n, double min, double max);


int main(int argc, char** argv){
    ifstream fileIn("dati.dat");
    ofstream fileOut("risultati.dat");
    double f, F, a, d;
    double stdev_f, stdev_F, stdev_a, stdev_d;
    double sigma_sist;
    result* vc = new result[N_SETS];
    set measures[N_SETS];
    double** dif_omega = new double*[N_SETS];
    double** dif_delta = new double*[N_SETS];
    double** c_values = new double*[N_SETS];
    int n[N_SETS];

    if(argc<N_SETS+1){
        cout << "Uso del programma: ./laboefm_c <n_misure_orario>  ";
        cout << "<n_misure_antiorario>  <n_misure_max_escursione>" << endl << endl;
        return -1;
    }

    srand(time(NULL));

    //Lettura
    ignoreText(&fileIn);
    fileIn >> f;
    ignoreText(&fileIn);
    fileIn >> stdev_f;
    ignoreText(&fileIn);
    fileIn >> F;
    ignoreText(&fileIn);
    fileIn >> stdev_F;
    ignoreText(&fileIn);
    fileIn >> a;
    ignoreText(&fileIn);
    fileIn >> stdev_a;
    ignoreText(&fileIn);
    fileIn >> d;
    ignoreText(&fileIn);
    fileIn >> stdev_d;

    for(int i=0; i<N_SETS; i++){
        n[i] = atoi(argv[i+1]);
    
        measures[i].dim = n[i];
        measures[i].v = new measure[n[i]];
        if(i==0 or i==1) n[i]--;
        else if(i==2) n[i] = n[i]/2;
        dif_delta[i] = new double[n[i]];
        dif_omega[i] = new double[n[i]];
        ignoreText(&fileIn);

        for(int k=0; k<measures[i].dim; k++){
            readMeasure(&fileIn, measures[i].v + k);
        }

        if(i==0 or i==1) reldif(measures[i], dif_omega[i], dif_delta[i], n[i]);
        else if(i==2) reldif2(measures[i], dif_omega[i], dif_delta[i], n[i]);
    }

    fileIn.close();

    //Stampa
    fileOut << "f = " << f << " +- " << stdev_f << " m" << endl;
    fileOut << "F = " << F << " +- " << stdev_F << " m" << endl;
    fileOut << "a = " << a << " +- " << stdev_a << " m" << endl;
    fileOut << "d = " << d <<  " +- " << stdev_d << " m" << endl << endl;
    fileOut << "Misure: omega (rad/s), delta (m)" << endl;
    
    for(int i=0; i<N_SETS; i++){
        printSet(&fileOut, measures[i]);
        fileOut << endl;

        for(int k=0; k<n[i]; k++){
            fileOut << fixed << showpoint << setprecision(2) << setw(8) << dif_omega[i][k];
            fileOut << fixed << showpoint << setprecision(5) << setw(10) << dif_delta[i][k] << endl;
        }
        fileOut << '(' << argv[i+1] << ')' << endl << endl;
    }

    //Medie e devst
    fileOut << scientific;

    for(int i=0; i<N_SETS; i++){
        c_values[i] = new double[n[i]];

        fileOut << "Valori di c, set n. " << i+1 << ":" << endl;

        for(int k=0; k<n[i]; k++){
            c_values[i][k] = 4*F*d*d*dif_omega[i][k] / ((d+a-F)*dif_delta[i][k]);
            fileOut << setprecision(3) << c_values[i][k] << endl;
        }

        setResult(vc[i], c_values[i], n[i]);

        fileOut << endl << "Media e deviazione standard: " << setprecision(3) << vc[i].m << " +- " << vc[i].stdev << endl << endl;
    }


    result c(vc, N_SETS);

    fileOut << "Deviazione standard della media (pesata): " << c.stdev << endl;

   
    //Propagazione degli errori 
    sigma_sist = sqrt(pow(stdev_a*dc_da(c.m, F, a, d), 2) + pow(stdev_d*dc_dd(c.m, F, a, d), 2) + pow(stdev_F*dc_dF(c.m, F, a, d) , 2));

    fileOut << setprecision(3) << "Errore \"sistematico\": " << sigma_sist << " m/s" << endl << endl;

    //Risultato
    c.stdev = sqrt(pow(sigma_sist, 2) + pow(c.stdev, 2));

    fileOut << "Misura di c:" << endl;
    fileOut << setprecision(3) << "c = " << c.m << " +- " << setprecision(1) << c.stdev << " m/s" << endl;

    double z = abs(c.m - C)/c.stdev;
    double p = 1 - stGauusianIntegral(N_P, -z, z);
    if(p<0.01) fileOut << "p < 1%" << endl << endl;
    else fileOut << "p = " << fixed << showpoint << setprecision(1) << p*100 << "%" << endl;

    //Termine programma
    fileOut.close();
    cout << "Analisi dati in \"risultati.dat\"." << endl;

    for(int i=0; i<N_SETS; i++){
        delete[] dif_omega[i];
        delete[] dif_delta[i];
        delete[] c_values[i];
        delete[] measures[i].v;
    }
    delete[] dif_omega;
    delete[] dif_delta;
    delete[] c_values;
    delete[] vc;

    return 0;
}

void ignoreText(ifstream* is){
    int c;

    *is >> ws;
    c = is->peek();
    while(!isdigit(c) and c!='-'){
        is->seekg(1,ios_base::cur);
        *is >> ws;
        c = is->peek();
    }
}

void readMeasure(ifstream* is, measure* m){
    double freq, dist;
    *is >> freq;
    m->nu = freq;
    m->omega = 2*M_PI*m->nu;
    *is >> dist;
    m->delta = dist/1000;
}

void printMeasure(ofstream* os, const measure m){
    *os << fixed << showpoint << setprecision(2) << setw(8) << m.omega;
    *os << fixed << showpoint << setprecision(5) << setw(10) << m.delta << endl;
}

void printSet(ofstream *os, const set s){
    for(int k=0; k<s.dim; k++){
        printMeasure(os, s.v[k]);
    }
}

void reldif(const set s, double* omega, double* delta, int n){
    for(int k=0; k<n; k++){
        omega[k] = s.v[k+1].omega - s.v->omega;
        delta[k] = s.v[k+1].delta - s.v->delta;
    }
}

void reldif2(const set s, double* omega, double* delta, int n){
    for(int k=0; k<n; k++){
        omega[k] = s.v[2*k].omega - s.v[2*k+1].omega;
        delta[k] = s.v[2*k].delta - s.v[2*k+1].delta;
    }
}



result::result(){
    m = 0;
    stdev = 0;
}

result::result(double* v, int n){
    double sum = 0;
    double sum2 = 0;

    for(int i=0; i<n; i++){
        sum += v[i];
        sum2 += v[i]*v[i];
    }

    m = sum/n;
    stdev = sqrt(sum2/n - m*m)*sqrt(n/(n-1));
}

result::result(result* vr, int n){
    double* w = new double[n];
    double sumpx = 0;
    double sump = 0;

    for(int i=0; i<n; i++){
        w[i] = pow(vr[i].stdev, -2);
        sumpx += w[i]*vr[i].m;
        sump += w[i];
    }

    m = sumpx/sump;
    stdev = pow(sump, -0.5);

    delete[] w;
}

void setResult(result& r, double* v, int n){
    double sum = 0;
    double sum2 = 0;

    for(int i=0; i<n; i++){
        sum += v[i];
        sum2 += v[i]*v[i];
    }

    r.m = sum/n;
    r.stdev = sqrt(sum2/n - r.m*r.m)*sqrt(n/(n-1));
}



double dc_da(double c, double F, double a, double d){
    double domega_ddelta_best;
    double result;

    domega_ddelta_best = domega_ddelta(c, F, a, d);   
    result = -4*domega_ddelta_best * F*d*d / pow((d+a-F), 2);

    return result;
}

double dc_dd(double c, double F, double a, double d){
    double domega_ddelta_best;
    double result;

    domega_ddelta_best = domega_ddelta(c, F, a, d); 
    result = 4*domega_ddelta_best * F*d*(d + 2*a -2*F) / pow((d+a-F), 2);

    return result;
}

double dc_dF(double c, double F, double a, double d){
    double domega_ddelta_best;
    double result;

    domega_ddelta_best = domega_ddelta(c, F, a, d); 
    result = 4*domega_ddelta_best * d*d*(d+a) / pow((d+a-F), 2);

    return result;
}

double domega_ddelta(double c, double F, double a, double d){
    double result;

    result = c * (d+a-F) / (4*F*d*d);

    return result;
}



double rand_double(double min, double max){
    double num;
    double x;

    x = (double)rand()/RAND_MAX;
    num = min + x*(max-min);

    return num;
}

double stGaussian(double z){
    double result;
    result = 1./sqrt(2*M_PI)*pow(M_E, -(z*z)/2);

    return result;
}

double stGauusianIntegral(int n, double min, double max){
    double sum = 0;
    double x;

    for (int i = 0; i < n; i++){
        x=rand_double(min, max);
        sum += stGaussian(x);
    }
    return sum*abs(max-min)/n;
}
