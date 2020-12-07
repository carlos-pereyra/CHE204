#include <math.h>
#include <iostream>
#include "poisson2d.h"

#ifndef DBG
#define DBG 1
#endif

using namespace std;

Poisson2D::Poisson2D(int n, int m) {
    if (DBG) printf("Poisson2D::Poisson2D()\n");
    Uold = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    U = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    x = (double *) malloc(sizeof(double)*n);    //real(n) vector
    y = (double *) malloc(sizeof(double)*m);    //real(m) vector
    f = (double *) malloc(sizeof(double)*n*m);  //real(nxm) vector
    mem=0;
    mem += sizeof(double)*n*m;
    mem += sizeof(double)*n*m;
    mem += sizeof(double)*n;
    mem += sizeof(double)*m;
    mem += sizeof(double)*n*m;

    nelem=n;
    melem=m;
    init();
}


Poisson2D::~Poisson2D() {
    free(Uold);
    free(U);
    free(x);
    free(y);
    free(f);
}

void Poisson2D::init() {
    if (DBG) printf("Poisson2D::init()\n");
    if (DBG) printf("n = %d\n", nelem);
    if (DBG) printf("m = %d\n", melem);
    // fill x with positions, domain omega(0,1)

    // file y with positions, domain omega(0,1)


    // fill U with intial guess
    for(int j=0; j<nelem; j++) {
        for(int i=0; i<melem; i++) {
            U[i + j*nelem] = 1;
        }
    }
    // fill f with initial boundary conditions
    for(int j=0; j<melem; j++) {
        for(int i=0; i<nelem; i++) {
            if(i==0 || i==(nelem-1)) { f[i + j*nelem] = 1; } 
            else { f[i + j*nelem] = 0; }
        }
    }
}

void Poisson2D::smooth(double *v) {
    /* jacobi relaxation smoothy algorithm.
    *
    *                up
    *                |
    *                |
    *   left ----- center ----- right
    *                |
    *                |
    *                down
    *
    */
    double up, down, right, left;
    int n=nelem;
    int m=melem;

    for(int j=1; j<(m-1); j++) {
        for(int i=1; i<(n-1); i++) {
            // save old field
            Uold[i + j*n] = U[i + j*n]
        }
    }

    for(int j=1; j<(m-1); j++) {
        for(int i=1; i<(n-1); i++) {
            // laplacian finite difference
            up          =U[i + (j+1)*n];
            down        =U[i + (j-1)*n];
            right       =U[(i+1) + j*n];
            left        =U[(i-1) + j*n];
            v[i + j*n]  =f[i + j*n] - (up + down + right + left) / 4;
            v[i + j*n]  /= 4;
            U[i + j*n]  = v[i + j*n];
        }
    }

}

void Poisson2D::writematrix2file(int index, double* matrix, string filename) {
    /* 
    assumes matrix is real(nxm) 
    */
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
