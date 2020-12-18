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
//      usage:
//          make
//          ./poisson
//
//      files:
//          mkdir data

#define NELEM 64
#define MELEM 64
#define ITERATIONS 2000
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
    double* u   = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    double* p   = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    
    double* dl  = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    double* rdl = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    
    double* utmp= (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector

    // --------- 2D POISSON SOLVER OPERATIONS ---------
    Poisson2D* ps2d = new Poisson2D(n,m, u,p);
    int l=log2(n);
    int nu=0;
    
    ps2d->copy(l, u,utmp); // duplicate u

    for(int index=0; index<v; index++) {
        // jacobi method
        ps2d->smooth(l,nu, u,p,"relax");
        ps2d->smooth(l,nu, utmp,p,"cgc");

        // coarse-grid correction
        ps2d->defect(l,nu, utmp,p,dl);          // risidual on l'th level
        ps2d->restriction(l-1,nu+1, dl,rdl);    // l'th level -> l'th-1 level
        ps2d->prolongation(l,nu, rdl,dl,utmp);  // l'th-1 level -> l'th level
        
        ps2d->smooth(l,nu, utmp,p,"cgc");
        
        // output
        std::string filename_utmp   = "data/pottmp" + to_string(index) + ".dat";
        std::string filename_u      = "data/pot" + to_string(index) + ".dat";
        std::string filename_errtmp = "data/errtmp" + to_string(index) + ".dat";
        std::string filename_err    = "data/err" + to_string(index) + ".dat";
        std::string filename_p      = "data/charge" + to_string(index) + ".dat";
        ps2d->index = index;
        ps2d->writematrix2file(filename_utmp,  "potentialtmp");
        ps2d->writematrix2file(filename_u,     "potential");
        ps2d->writematrix2file(filename_errtmp,"errortmp");
        ps2d->writematrix2file(filename_err,   "error");
        ps2d->writematrix2file(filename_p,     "charge");

        // error residual
        printf("l=%d nu=%d iteration=%d error= %lf errortmp=%lf\n",l, nu, index, ps2d->error, ps2d->errortmp);
        //if(ps2d->errortmp<1e-4) {
        if(ps2d->error<1e-4) {
            break;
        }
    }

    free(u);
    free(p);
    free(dl);
    free(rdl);
    free(utmp);
    delete ps2d;
}
