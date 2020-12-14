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
// library
#include "poisson2d.h"

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
//      usage:
//          make
//          ./poisson
//

#define NELEM 51
#define MELEM 51
#define ITERATIONS 2000
using namespace std;

int main(int argc, char** argv) {
    int n=NELEM, m=MELEM;
    int v=ITERATIONS;
    
    // --------- 2D POISSON SOLVER OPERATIONS ---------
    //printf("start program\n");
    //printf("n = %d\n", n);
    //printf("m = %d\n", m);
    //printf("v = %d\n", v);
    Poisson2D* ps2d = new Poisson2D(n,m);

    // --------- POTENTIAL VECTOR ---------
    double* u = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    double* p = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector

    for(int index=0; index<v; index++) {
        // jacobi method
        ps2d->smooth(u);

        // multigrids
        //ps2d->residual();

        // error residual
        printf("iteration= %d error= %lf\n", index, ps2d->error);
        if(ps2d->error<1e-4) {
            break;
        }

        // output
        std::string filenameu = "data/potential" + to_string(index) + ".dat";
        std::string filenamep = "data/charge" + to_string(index) + ".dat";
        ps2d->writematrix2file(filenameu,"potential");
        ps2d->writematrix2file(filenamep,"charge");
    }

    free(u);
    free(p);
    delete ps2d;
}
