#include <cstdlib>    // c++
#include <cstdio>     // c++
#include <omp.h>
#include <vector>     // c++
#include <algorithm>  // c++
#include <iostream>   // c++
#include <iomanip>    // c++
#include <fstream>    // c++
#include <math.h>
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.c
//
//      ./test
//

#define ranfloat(w) ( w * ((rand() / (double) RAND_MAX)  - 0.5) )
#define DBG 0
using namespace std;

float f(double x, int n, int m, double *rx);
float series(long trials);
void write2file(const long n, float result, float mcresult);

int main(int argc, char** argv) {
    double x, y, z;
    double rx[2] = {2.4, 1.2};
    double ry[2] = {0, 2.078};


    // parameters
    if(DBG) printf("M_PI = %f\n", M_PI);
    // fourier series sum 
    for(n=-4; n <= 4; n++) {
        for(m=-4; m <= 4; m++) {
            z += cos(kx*n*x)
        }
    }
    write2file(x, y, z);
}

double f(double x, int n, int m, double* rx) {
    double sigma = 0.1;
    return exp( -(x-n*rx[0]-m*rx[1]) / sigma );
}

float montecarlo_integration(long trials) {
    float ratio = 1;
    return ratio;
}

void write2file(const long n, float result, float mcresult) {
    string filename = "data/output.dat";
    ofstream outfile;
    if (n==100) {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# n";
        outfile << left << setw(12) << "result";
        outfile << left << "mc-result\n";
    } else {
        outfile.open(filename, std::fstream::app);
    }
    // print values
    outfile << std::scientific;
    outfile.precision(4);
    outfile << left << setw(12) << n;
    outfile << left << setw(12) << result;
    outfile << left << mcresult << "\n";

    outfile.close();
}

