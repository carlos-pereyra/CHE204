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

float f(float x, float y);
float montecarlo_integration(long trials);
void write2file(const long n, float result, float mcresult);

int main(int argc, char** argv) {
    long p, n, i; long nstart = 100;
    float x, y, z, s, mcresult;
    float result=(pow(M_PI, 2) / 2);
    // parameters
    if(DBG) printf("RAND_MAX = %d\n", RAND_MAX);
    if(DBG) printf("M_PI = %f\n", M_PI);
    if(DBG) printf("upper limit = %f\n", M_PI/2);
    if(DBG) printf("lower limit = %f\n\n", -M_PI/2);
    // p'th order (mc steps loop)
    for(p=1; p < 8; p++) {
        n = pow(nstart, p/float(2));
        s = 0;
        // monte carlo integration loop
        for(i=0; i<n; i++) {
            x = ranfloat(M_PI);
            y = ranfloat(M_PI);
            // sum
            s += f(x,y);
        }
        mcresult = (s / n) * M_PI * M_PI; // box area = pi^2
        printf("n = %ld\n", n);
        printf("mc result = %f\n", mcresult);
        printf("analytical result = %f\n", result);
        write2file(n, result, mcresult);
    }
}

float f(float x, float y) {
    return pow(cos(3*x), 2) + 2 * y * pow(sin(1.3*y), 2);
}

float montecarlo_integration(long trials) {
    float ratio = 1;
    return ratio;
}

void write2file(const long n, float result, float mcresult) {
    string filename = "data/output.dat";
    ofstream outfile;
    outfile.open(filename, std::fstream::out);
    outfile << left << setw(12) << "# n";
    outfile << left << setw(12) << "result";
    outfile << left << "mc-result\n";
    // print values
    outfile << std::scientific;
    outfile.precision(4);
    outfile << left << setw(12) << n;
    outfile << left << setw(12) << result;
    outfile << left << mcresult;

    outfile.close();
}
