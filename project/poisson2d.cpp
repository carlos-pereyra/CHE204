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

#define SCALEX 1
#define SCALEY 1

using namespace std;

Poisson2D::Poisson2D(int n, int m) {
    if (DBG) printf("Poisson2D::Poisson2D()\n");
    uold = (double *) malloc(sizeof(double)*n*m);   //real(nxm) vector
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
    init();
}


Poisson2D::~Poisson2D() {
    free(uold);
    free(u);
    free(x);
    free(y);
    free(p);
}

void Poisson2D::init() {
    if (DBG) printf("\nPoisson2D::init()\n");
    if (DBG) printf(" n = %d\n", nelem);
    if (DBG) printf(" m = %d\n", melem);
    
    // fill u (electrostatic-potential) with intial guess
    for(int j=0; j<nelem; j++) {
        for(int i=0; i<melem; i++) {
            u[i + j*nelem] = 0;
        }
    }

    // fill p (charge density) with initial boundary conditions
    for(int j=0; j<melem; j++) {
        for(int i=0; i<nelem; i++) {
            //if(i==1 || i==(nelem-2)) { p[i + j*nelem] = 1; } 
            if(i==ceil((nelem-2)*0.5)) { p[i + j*nelem] = 1; } 
            else { p[i + j*nelem] = 0; }
        }
    }

    // finite difference variables
    dx   = SCALEX * 1 / (double) nelem;
    dy   = SCALEY * 1 / (double) melem;
    dx2  = pow(dx, 2);
    dy2  = pow(dy, 2);

}

double Poisson2D::smooth(double *v) {
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
    *
    *   return percent difference betwee k'th and k+1'th iteration
    *
    */
    double up, down, right, left;
    int n=nelem;
    int m=melem;

    double trace = 0; double traceold = 0;
    for(int j=1; j<(m-1); j++) {
        for(int i=1; i<(n-1); i++) {
            // laplacian finite difference
            up          =u[i + (j+1)*n];
            down        =u[i + (j-1)*n];
            right       =u[(i+1) + j*n];
            left        =u[(i-1) + j*n];
            v[i + j*n]  =0.5 * ( dx2*(up + down)/(dx2 + dy2) + dy2*(right + left)/(dx2 + dy2) + (dx2*dy2)*p[i + j*n]/(dx2 + dy2) );
            // error sum
            if(i==j) {
                trace    +=v[i + j*n];
                traceold +=u[i + j*n];
            }
            // update private potential u
            u[i + j*n]  = v[i + j*n];

        }
    }

    error = abs(pow(trace,2) - pow(traceold,2)) / pow(traceold,2);
    return error;
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
    
    //if (i==0) {
    outfile.open(filename, std::fstream::out);
    outfile << left << setw(12) << "# x";
    outfile << left << setw(12) << "y";
    outfile << left << "f(x,y)" << "\n";
    //else if(newlineflag) {
    //    outfile.open(filename, std::fstream::app);
    //    outfile << "\n";
    //    outfile.close();}
    //else {
    //    outfile.open(filename, std::fstream::app);
    //}
    
    // print values
    if(mode=="potential") {
        outfile << std::scientific;
        outfile.precision(4);
        for(int j=1; j<(m-1); j++) {
            for(int i=1; i<(n-1); i++) {
                outfile << left << setw(12) << i;
                outfile << left << setw(12) << j;
                outfile << left << u[i + j*n] << "\n";
            }
            outfile << "\n";
        }
    }
    else if (mode=="charge") {
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
