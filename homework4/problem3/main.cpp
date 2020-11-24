/*
 * Made by Carlos Pereyra Ph.D. student, UC Davis, CHE 204 2020
 */
#include <cstdlib>    // c++
#include <cstdio>     // c++
//#include <omp.h>
#include <vector>     // c++
#include <algorithm>  // c++
#include <iostream>   // c++
#include <iomanip>    // c++
#include <fstream>    // c++
#include <math.h>
#include <sstream>
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.cpp
//
//      ./test
//

//#define N 60
//#define NX 60
//#define NY 60
//#define SCALE 10
#define DUTY 0.9
using namespace std;

double** inputsignal(string filename, double **addr, double *xj, double *yj, int *n);
double w(int j, int n);
double gauss(int j, int n);
void write2file(string filename, int i, int newlineflag, double x, double y, double z, double q);

int main(int argc, char** argv) {
    double freal, fimag, fabs;
    double inv_Freal, inv_Fimag, inv_Fabs;
    double inv_Freal_w, inv_Fimag_w, inv_Fabs_w;
    double inv_Freal_gauss, inv_Fimag_gauss, inv_Fabs_gauss;

    double **address = (double**) malloc(sizeof(double*)*2); // points to xj and yj
    double *xj = (double*) malloc(sizeof(double)*1);
    double *yj = (double*) malloc(sizeof(double)*1);
    
    // input signal
    int n;
    inputsignal("noisyspectrum.dat", address, xj, yj, &n);
    xj = address[0];
    yj = address[1];

    // parameters
    printf("PARAMETERS\n");
    printf("===============\n");
    printf("M_PI = %lf\n", M_PI);
    printf("n = %d\n", n);
    printf("midpoint n/2 = %d\n", n/2);

    // nk = nj = n
    double *ck_real = (double*) malloc(sizeof(double)*n);
    double *ck_imag = (double*) malloc(sizeof(double)*n);
    double theta   = 2*M_PI / n;

    // fourier transform
    for(int k=0; k<n; k++) {
        freal = 0;
        fimag = 0;
        //inv_Freal = 0;
        //inv_Fimag = 0;
        for(int j=0; j<n; j++) {
            // FT
            freal += yj[j]*cos(k*theta*j);  // + yjimag[j]*sin(k*theta*j); // (if there were an imaginary part)
            fimag += -yj[j]*sin(k*theta*j); // + yjimag[j]*sin(k*theta*xj[j]); // (if there were an imaginary part)

            // inv FT
            //inv_Freal += (cos(xj[k]*theta*j)*(yj[k]*cos(xj[k]*theta*j)) + (yj[k]*sin(xj[k]*theta*j))*sin(xj[k]*theta*j))* w(j,n);
            //inv_Fimag += (sin(xj[k]*theta*j)*(yj[k]*cos(xj[k]*theta*j)) - (yj[k]*sin(xj[k]*theta*j))*cos(xj[k]*theta*j))* w(j,n);
        }
        // FT
        fabs = sqrt(pow(freal, 2) + pow(fimag, 2));
        ck_real[k] = freal / sqrt(2*M_PI);
        ck_imag[k] = fimag / sqrt(2*M_PI);
        // inverse FT
        //inv_Fabs = sqrt(pow(inv_Freal, 2) + pow(inv_Fimag, 2));
        // write [filename] [index] [newlineflag] [val1] [val2] [val3] [val4]
        write2file("data/signal.dat", k, 0, xj[k], yj[k], 0, 0);
        write2file("data/ft.dat", k, 0, k, ck_real[k], k, ck_imag[k]);
        //write2file("data/filteredsignal.dat", k, 0, xj[k], fk[k], 0, 0);
    }

    // inverse fourier transform
    for(int k=0; k<n; k++) {
        inv_Freal       = 0;    // original - inverse fourier transform (real component)
        inv_Fimag       = 0;    // original - inverse fourier transform (imaginary component)
        inv_Freal_w     = 0;    // heaviside filtered - inverse fourier transform (real component)
        inv_Fimag_w     = 0;    // heaviside filtered - inverse fourier transform (imaginary component)
        inv_Freal_gauss = 0;    // gaussian filtered - inverse fourier transform (real component)
        inv_Fimag_gauss = 0;    // gaussian filtered - inverse fourier transform (imaginary component)
        //for(int j=(-n/2); j<(n/2); j++) {
        for(int j=0; j<n; j++) {
            double phase = k*theta*j;
            double Re = (cos(phase)*ck_real[j] - sin(phase)*ck_imag[j]);
            double Im = (sin(phase)*ck_real[j] + cos(phase)*ck_imag[j]);
            inv_Freal       += Re;
            inv_Fimag       += Im; 
            inv_Freal_w     += Re * w(j,n); 
            inv_Fimag_w     += Im * w(j,n);
            inv_Freal_gauss += Re * gauss(j,n); 
            inv_Fimag_gauss += Im * gauss(j,n);
        }
        // magnitude
        inv_Fabs        = sqrt(pow(inv_Freal, 2) + pow(inv_Fimag, 2));
        inv_Fabs_w      = sqrt(pow(inv_Freal_w, 2) + pow(inv_Fimag_w, 2));
        inv_Fabs_gauss  = sqrt(pow(inv_Freal_gauss, 2) + pow(inv_Fimag_gauss, 2));
        // normalize
        inv_Fabs        *= sqrt(2*M_PI) / n;
        inv_Freal       *= sqrt(2*M_PI) / n;
        inv_Fimag       *= sqrt(2*M_PI) / n;

        inv_Fabs_w      *= sqrt(2*M_PI) / (DUTY * n);
        inv_Freal_w     *= sqrt(2*M_PI) / (DUTY * n);
        inv_Fimag_w     *= sqrt(2*M_PI) / (DUTY * n);

        inv_Fabs_gauss  *= sqrt(2*M_PI) / n;
        inv_Freal_gauss *= sqrt(2*M_PI) / n;
        inv_Fimag_gauss *= sqrt(2*M_PI) / n;
        // write [filename] [index] [newlineflag] [val1] [val2] [val3] [val4]
        write2file("data/inv_ft_abs.dat", k, 0, xj[k], inv_Fabs, 0, 0);
        write2file("data/inv_ft_w_abs.dat", k, 0, xj[k], inv_Fabs_w, 0, 0);
        write2file("data/inv_ft_gauss_abs.dat", k, 0, xj[k], inv_Fabs_gauss, 0, 0);
        //
        write2file("data/inv_ft.dat", k, 0, xj[k], inv_Freal, xj[k], inv_Fimag);
        write2file("data/inv_ft_w.dat", k, 0, xj[k], inv_Freal_w, xj[k], inv_Fimag_w);
        write2file("data/inv_ft_gauss.dat", k, 0, xj[k], inv_Freal_gauss, xj[k], inv_Fimag_gauss);
    }

    //double complex pk[n];
    //int LENSAVE = 2*n + int(log(n)/log(2)) + 4;
    //double wsave[n];
    //int LENWORK = 2*n;
    //double work[n];
    //_cfft1b(n, 1, &pk, n, &wsave, LENSAVE, &work, );

    // free memory
    free(address[0]);
    free(address[1]);
    free(address);
    free(ck_real);
    free(ck_imag);
}

// INPUT

double** inputsignal(string filename, double **addr, double *xj, double *yj, int *n) {
    // read in positions from .xyz file
    ifstream infile(filename);
    string line;
    double x, y;
    int index = 0;

    // IO
    if(!infile.is_open()) {
        printf("input file not found\n");
        exit(1);
    }
    //std::getline(infile, line); std::istringstream iss(line);
    //std::cout << "natoms = " << natoms << "\n";
    //std::getline(infile, line);  // line 2

    // READ LINES
    while(std::getline(infile, line)) {
        // read x y z coords
        std::istringstream iss(line);

        if (iss >> x >> y) //cout << "line = "<< index << " x = " << x << " y = " << y << "\n";

        xj[index] = double(x);
        yj[index] = double(y);
        
        index++;
        double* tmp_xj = (double* ) realloc(xj, (index+1)*sizeof(double));
        double* tmp_yj = (double* ) realloc(yj, (index+1)*sizeof(double));
        if(tmp_xj != NULL) xj = tmp_xj;
        if(tmp_yj != NULL) yj = tmp_yj;
    }
    infile.close();
    
    addr[0] = xj;
    addr[1] = yj;
    n[0] = index;

    return addr;
}

// WINDOW FUNCTION

double w(int j, int n) {
    double value;
    if(j <= int(DUTY * n)) value = 1;
    else value = 0;
    return value;
}

double gauss(int j, int n) {
    double sigma = 0.1;
    return exp(-0.5*pow(2*n / (sigma * (n-1)), 2));
}

// OUTPUT

void write2file(string filename, int i, int newlineflag, double x, double y, double z, double q) {
    ofstream outfile;
    if (i==0) {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# x";
        outfile << left << setw(12) << "y";
        outfile << left << setw(12) << "z";
        outfile << left << "q\n";
    } else if(newlineflag) {
        outfile.open(filename, std::fstream::app);
        outfile << "\n";
        outfile.close();
    } else {
        outfile.open(filename, std::fstream::app);
    }
    // print values
    outfile << std::scientific;
    outfile.precision(4);
    outfile << left << setw(12) << x;
    outfile << left << setw(12) << y;
    outfile << left << setw(12) << z;
    outfile << left << q << "\n";
    outfile.close();
}
