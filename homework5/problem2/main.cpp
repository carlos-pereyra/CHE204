#include <cstdlib>    // c++
#include <cstdio>     // c++
#include <omp.h>
#include <vector>     // c++
#include <algorithm>  // c++
#include <iostream>   // c++
#include <iomanip>    // c++
#include <fstream>    // c++
#include <math.h>
#include <mpi.h>
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.cpp
//
//      ./test
//
//      for extra juiciness use MPI - currently does not work :(
//      
//      mpicxx main.cpp -o test
//
//      mpirun -n 2 ./test
//

#define N 1000
#define DT 0.01
using namespace std;

double F1(double u, double v);
double F2(double u, double v);
void write2file(string filename, int i, double x, double y, double z);

int main(int argc, char** argv) {
    double uk = 1, ukHalf, ukNew; // initial condition u0 = 1
    double vk = 1, vkHalf, vkNew; // initial condition v0 = 1
    double f1, f2;
    double errU, errV;
    double t; 

    int n = N;
    double dt = DT;

    for(int i=0; i<n; i++) {
        t = dt * i;
        // output
        write2file("data/output.dat", i, t, uk, vk);
        printf("\nIteration=%d\n", i);
        printf("U=%lf V(t)=%lf\n", uk, vk);
        // half step calculation
        ukHalf = uk + 0.5 * dt * F1(uk,vk);
        vkHalf = vk + 0.5 * dt * F2(uk,vk);
        // new u, new v
        ukNew = uk + dt * F1(ukHalf, vkHalf);
        vkNew = vk + dt * F2(ukHalf, vkHalf);
        // error
        errU = abs(uk - ukNew);
        errV = abs(vk - vkNew);
        uk = ukNew;
        vk = vkNew;
    }

    printf("done! \n");
}

//
// SYSTEM OF ODE'S
//

double F1(double u, double v) {
    // function
    //
    //      F1(u,v) = du/dt = (a*u - b*u**2) - c*u*v
    //
    double a = 4;
    double b = 2;
    double c = 1;
    double val = (a*u - b*pow(u,2)) - c*u*v;
    return val;
}

double F2(double u, double v) {
    // function
    //
    //      F2(u,v) = dv/dt = -d*v + e*u*v
    //
    double d = 2;
    double e = 2;
    double val = -d*v + e*u*v;
    return val;
}

//
// OUTPUT
//

void write2file(string filename, int i, double x, double y, double z) {
    ofstream outfile;
    if (i==0) {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# t";
        outfile << left << setw(12) << "u(t)";
        outfile << left << setw(12) << "v(t)";
    }
    //else if(newlineflag) {
    //    outfile.open(filename, std::fstream::app);
    //    outfile << "\n";
    //    outfile.close();
    //} 
    else {
        outfile.open(filename, std::fstream::app);
    }
    // print values
    outfile << std::scientific;
    outfile.precision(4);
    outfile << left << setw(12) << x;
    outfile << left << setw(12) << y;
    outfile << left << z << "\n";
    outfile.close();
}
