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

#define NELEM 64 // both have to be base 2
#define MELEM 64
#define ITERATIONS 500
using namespace std;

int main(int argc, char** argv) {
    int n=NELEM, m=MELEM;
    int v=ITERATIONS;
    
    /* --------- VECTOR INITIALIZATION ---------
     *
     * poisson:
     *   \nabla^2 u = p
     *
     */
    double* u  = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    double* p  = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    double* dl = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    double* tmp= (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector

    // --------- 2D POISSON SOLVER OPERATIONS ---------
    Poisson2D* ps2d = new Poisson2D(n,m, u,p);

    int l=log2(n);
    int nu=0;
    for(int index=0; index<v; index++) {
        // jacobi method
        ps2d->smooth(l,nu, u,p);

        // coarse-grid correction
        //ps2d->defect(l,nu, u,p,dl);
        //ps2d->restriction();
        //ps2d->prolongation();

        // output
        std::string filenameu = "data/potential" + to_string(index) + ".dat";
        std::string filenamep = "data/charge" + to_string(index) + ".dat";
        std::string filenamee = "data/error" + to_string(index) + ".dat";
        ps2d->writematrix2file(filenameu,"potential");
        ps2d->writematrix2file(filenamep,"charge");
        ps2d->write2file(filenamee,index,ps2d->error,0);

        // error residual
        printf("iteration= %d error= %lf\n", index, ps2d->error);
        if(ps2d->error<1e-3) {
            break;
        }
    }

    free(u);
    free(p);
    free(dl);
    free(tmp);
    delete ps2d;
}
