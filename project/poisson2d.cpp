#include <math.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>      // file input, output
#include <iomanip>      // print setprecision 
#include "poisson2d.h"

#ifndef DBG
#define DBG 0
#endif

#define LENGTHB 1       // BAR LENGTH
#define LENGTHX 3       // X DIRECTION LENGTH
#define LENGTHY 3       // Y DIRECTION LENGTH

using namespace std;

Poisson2D::Poisson2D(int n, int m, double* phi, double* rho) {
    if (DBG) printf("Poisson2D::Poisson2D()\n");
    u = (double *) malloc(sizeof(double)*n*m);      //real(nxm) vector
    x = (double *) malloc(sizeof(double)*n);        //real(n) vector
    y = (double *) malloc(sizeof(double)*m);        //real(m) vector
    p = (double *) malloc(sizeof(double)*n*m);      //real(nxm) vector
    mem=0;
    mem += sizeof(double)*n*m;
    mem += sizeof(double)*n*m;
    mem += sizeof(double)*n;
    mem += sizeof(double)*m;
    mem += sizeof(double)*n*m;

    nelem=n;
    melem=m;
    init(phi, rho);

    // update all phi and rho elements
    //phi = u;
    //rho = p;
}


Poisson2D::~Poisson2D() {
    free(u);
    free(x);
    free(y);
    free(p);
}

void Poisson2D::init(double* phi, double* rho) {
    if (DBG) printf("\nPoisson2D::init()\n");
    if (DBG) printf(" n = %d\n", nelem);
    if (DBG) printf(" m = %d\n", melem);
    
    // fill u (electrostatic-potential) with intial guess and 
    // dirichlet boundary conditions
    for(int j=0; j<melem; j++) {
        for(int i=0; i<nelem; i++) {
            if(DBG) printf("filling u[%d][%d]\n", i, j);
            if( (i==0) && (i==(nelem-1)) && (j==0) && (j==(nelem-1)) ) { 
                // boundary conditions
                phi[i + j*nelem] = u[i + j*nelem] = 0;
            } else {
                // internal guess
                phi[i + j*nelem] = u[i + j*nelem] = 0;
            }
        }
    }

    // fill p (charge density) with initial boundary conditions
    int N1 = ceil(0.5*melem*(1 - LENGTHB / (double) LENGTHY));
    int N2 = ceil(1.0*melem*(LENGTHB / (double) LENGTHY));
    for(int j=0; j<melem; j++) {
        for(int i=0; i<nelem; i++) {
            if(DBG) printf("filling rho[%d][%d]\n", i, j);

            // charged bar conditions
            if(i==ceil((nelem-2)*0.5) && j>=N1 && j<=(N1+N2)) { 
                rho[i + j*nelem] = p[i + j*nelem] = 1;
                //phi[i + j*nelem] = u[i + j*nelem] = 1;
            } 
            else { 
                rho[i + j*nelem] = p[i + j*nelem] = 0;
            }

        }
    }

    // finite difference
    dx   = LENGTHX * 1 / (double) nelem;
    dy   = LENGTHY * 1 / (double) melem;
    dx2  = pow(dx, 2);
    dy2  = pow(dy, 2);

}

void Poisson2D::smooth(double *v, double *f) {
    /* 
     * jacobi relaxation smoothy algorithm.
     *
     *                up
     *                |
     *                |
     *   left ----- center ----- right
     *                |
     *                |
     *                down
     *
     *
     * return percent difference.
     *
     */
    if (DBG) printf("\nPoisson2D::smooth()\n");
    
    double up, down, right, left;
    int n=nelem;
    int m=melem;

    double trace = 0;
    double traceold = 0;
    
    for(int j=1; j<(m-1); j++) {
        for(int i=1; i<(n-1); i++) {
            // laplacian finite difference
            up          =v[i + (j+1)*n];
            down        =v[i + (j-1)*n];
            right       =v[(i+1) + j*n];
            left        =v[(i-1) + j*n];
            u[i + j*n]  =v[i + j*n]; // u^n
            v[i + j*n]  =0.5 * ( dx2*(up + down)/(dx2 + dy2) +
                                 dy2*(right + left)/(dx2 + dy2) + 
                                 (dx2*dy2)*f[i + j*n]/(dx2 + dy2) ); // u^n+1
            // error
            if(i==j) {
                trace    +=v[i + j*n];
                traceold +=u[i + j*n];
            }

        }
    }

    error = (pow(trace,2) - pow(traceold,2)) / pow(traceold,2);
    if(DBG) printf("error = %lf\n", error);
}

void Poisson2D::defect(double *v, double *f, double* dl) {
    /*
     * defect calculation.
     *
     *      dl = Ll ul - fl
     *
     * determines +residual error, dl.
     */
    if (DBG) printf("\nPoisson2D::defect()\n");
    
    double up, down, right, left;
    int n=nelem;
    int m=melem;

    for(int j=1; j<(m-1); j++) { // interior defect points
        for(int i=1; i<(n-1); i++) {
            up          =v[i + (j+1)*n]; // Ll ul - fl operation
            down        =v[i + (j-1)*n];
            right       =v[(i+1) + j*n];
            left        =v[(i-1) + j*n];
            dl[i + j*n] =0.5 * ( dx2*(up + down)/(dx2 + dy2) +
                                 dy2*(right + left)/(dx2 + dy2) +
                                 (dx2*dy2)*f[i + j*n]/(dx2 + dy2) ) - f[i + j*n];

        }
    }

    for(int i=0; i<n; i++) dl[i + 0]   = dl[i + (m-1)*n] = 0; // boundary defect points
    for(int j=0; j<m; j++) dl[0 + j*n] = dl[(n-1) + j*n] = 0;

}

void Poisson2D::restriction(double *v, double *f) {
}

void Poisson2D::prolongation(double *v, double *f) {
}


void Poisson2D::writematrix2file(std::string filename, std::string mode) {
    /* 
    *
    * assumes matrix is real(nxm) 
    *
    */
    int newlineflag = 0;
    int n=nelem;
    int m=melem;
    ofstream outfile;
    
    if(mode=="potential") {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# x";
        outfile << left << setw(12) << "y";
        outfile << left << setw(12) << "f(x,y)";
        outfile << left << "error" << "\n";
        outfile << std::scientific;
        outfile.precision(4);
        for(int j=1; j<(m-1); j++) {
            for(int i=1; i<(n-1); i++) {
                outfile << left << setw(12) << i;
                outfile << left << setw(12) << j;
                outfile << left << setw(12) << u[i + j*n];
                outfile << left << error << "\n";
            }
            outfile << "\n";
        }
    }
    else if (mode=="charge") {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# x";
        outfile << left << setw(12) << "y";
        outfile << left << "f(x,y)" << "\n";
        outfile << std::scientific;
        outfile.precision(4);
        for(int j=0; j<m; j++) {
            for(int i=0; i<n; i++) {
                outfile << left << setw(12) << i;
                outfile << left << setw(12) << j;
                outfile << left << p[i + j*n] << "\n";
            }
            outfile << "\n";
        }
    }

}
