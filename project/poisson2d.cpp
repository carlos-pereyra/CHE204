#include <math.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>      // file input, output
#include <iomanip>      // print setprecision 
#include "poisson2d.h"

#ifndef DBG
#define DBG 1
#endif

#define LENGTHB 1       // BAR LENGTH
#define LENGTH 3        // SCALE SQUARE LATICE PLANE

using namespace std;

Poisson2D::Poisson2D(int n, int m, double* v, double* p) {
    if (DBG) printf("Poisson2D::Poisson2D()\n");
    u = (double *) malloc(sizeof(double)*n*m);      //real(nxm) vector
    x = (double *) malloc(sizeof(double)*n);        //real(n) vector
    y = (double *) malloc(sizeof(double)*m);        //real(m) vector
    f = (double *) malloc(sizeof(double)*n*m);      //real(nxm) vector
    mem=0;
    mem += sizeof(double)*n*m;
    mem += sizeof(double)*n*m;
    mem += sizeof(double)*n;
    mem += sizeof(double)*m;
    mem += sizeof(double)*n*m;

    nelem=n;
    melem=m;
    init(v, p);

}


Poisson2D::~Poisson2D() {
    free(u);
    free(x);
    free(y);
    free(f);
}

void Poisson2D::init(double* v, double* p) {
    if (DBG) printf("\nPoisson2D::init()\n");
    if (DBG) printf(" n = %d\n", nelem);
    if (DBG) printf(" m = %d\n", melem);
    
    // fill u (electrostatic-potential) with intial guess and 
    // dirichlet boundary conditions
    for(int j=0; j<melem; j++) {
        for(int i=0; i<nelem; i++) {
            //if(DBG) printf("filling u[%d][%d]\n", i, j);
            if( (i==0) && (i==(nelem-1)) && (j==0) && (j==(nelem-1)) ) { 
                // boundary conditions
                v[i + j*nelem] = u[i + j*nelem] = 0;
            } else {
                // internal guess
                v[i + j*nelem] = 0;
                u[i + j*nelem] = 0;
            }
        }
    }

    // fill p (charge density) with initial boundary conditions
    int N1 = ceil(0.5*melem*(1 - LENGTHB / (double) LENGTH));
    int N2 = ceil(1.0*melem*(LENGTHB / (double) LENGTH));
    for(int j=0; j<melem; j++) {
        for(int i=0; i<nelem; i++) {
            //if(DBG) printf("filling rho[%d][%d]\n", i, j);

            // charged bar conditions
            if(i==ceil((nelem-2)*0.5) && j>=N1 && j<=(N1+N2)) { 
                p[i + j*nelem] = f[i + j*nelem] = 1;
                //phi[i + j*nelem] = u[i + j*nelem] = 1;
            } 
            else { 
                p[i + j*nelem] = f[i + j*nelem] = 0;
            }

        }
    }

    // finite difference
    dx   = LENGTH * 1 / (double) nelem;
    dy   = LENGTH * 1 / (double) melem;
    dx2  = pow(dx, 2);
    dy2  = pow(dy, 2);

}

void Poisson2D::smooth(int l, int nu, double *v, double *p) {
    /* 
     * jacobi relaxation smoothy algorithm. iteratively solve \nabla^2 u = f
     * assumptions: (1) v, f are square lattices with 2^l+1 elements on each dimension
     *              (2) l is level number
     *              (3) nu is level depth number
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
    int n       =pow(2,l+1);
    int m       =pow(2,l+1);
    int step    =pow(2,nu);
    double d    =LENGTH / (double) n; 
    double d2   =pow(d, 2); //(double) pow(2,l+1);
    double d2i  =1 / d2;
    //dx2 = dy2 = d2;
    printf("step = %d nu = %d\n",step, nu);
    printf("d2i = %lf\n",d2i);

    double trace= 0;
    double traceold= 0;
    
    for(int j=step; j<(m-1); j+=step) {
        for(int i=step; i<(n-1); i+=step) {
            // laplacian finite difference
            up          =v[i + (j+1)*n];
            down        =v[i + (j-1)*n];
            right       =v[(i+1) + j*n];
            left        =v[(i-1) + j*n];
            u[i + j*n]  =v[i + j*n]; // u^n
            // general expression
            v[i + j*n]  =0.5 * ( dx2*(up + down)/(dx2 + dy2) +
                                 dy2*(right + left)/(dx2 + dy2) + 
                                 (dx2*dy2)*p[i + j*n]/(dx2 + dy2) ); // u^n+1
            // square expression 
            //v[i + j*n]  =0.25*((right + left + up + down) + d2i*f[i + j*n]);
            // error
            if(i==j) {
                trace    +=v[i + j*n];
                traceold +=u[i + j*n];
            }

        }
    }

    error = abs(pow(trace,2) - pow(traceold,2)) / pow(traceold,2);
    if(DBG) printf("error = %lf\n", error);
}

void Poisson2D::defect(int l, int nu, double *v, double *f, double* dl) {
    /*
     * defect calculation.
     *
     *      dl = Ll ul - fl
     *
     * determines +residual error, dl.
     */
    if (DBG) printf("\nPoisson2D::defect()\n");
    
    double up, down, right, left, onsite;
    int n       =pow(2,l+1);
    int m       =pow(2,l+1);
    int step    =pow(2,nu);
    double d2i  =1/(double) pow(2,l+1);

    for(int j=step; j<(m-1); j+=step) { // interior defect points
        for(int i=step; i<(n-1); i+=step) {
            up          =v[i + (j+1)*n]; // Ll ul - fl operation
            down        =v[i + (j-1)*n];
            right       =v[(i+1) + j*n];
            left        =v[(i-1) + j*n];
            onsite      =v[i + j*n];
            /*
            dl[i + j*n] =0.5 * ( dx2*(up + down)/(dx2 + dy2) +
                                 dy2*(right + left)/(dx2 + dy2) +
                                 (dx2*dy2)*f[i + j*n]/(dx2 + dy2) ) - f[i + j*n];*/

            // minus defect
            dl[i+j*n]   =-d2i*(right + left + up + down - 4*onsite) + f[i+j*n];
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
                outfile << left << f[i + j*n] << "\n";
            }
            outfile << "\n";
        }
    }

}

void Poisson2D::write2file(string filename, double x, double y, double z) {
    ofstream outfile;
    outfile.open(filename, std::fstream::out);
    outfile << left << setw(12) << "# x";
    outfile << left << setw(12) << "y";
    outfile << left << "z" << "\n";
    // print values
    outfile << std::scientific;
    outfile.precision(4);
    outfile << left << setw(12) << x;
    outfile << left << setw(12) << y;
    outfile << left << z << "\n";
    outfile.close();
}
