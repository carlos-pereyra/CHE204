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
#define ranfloat(w) ( w * ((rand() / (double) RAND_MAX)) ) // RANDOM NUMBER FROM 0 TO W

using namespace std;

Poisson2D::Poisson2D(int n, int m, double* phi, double* rho) {
    if (DBG) printf("Poisson2D::Poisson2D()\n");
    utmp = (double *) malloc(sizeof(double)*n*m);   //real(nxm) vector
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

}


Poisson2D::~Poisson2D() {
    free(utmp);
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

    // initial values
    
    // +charged bar conditions
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

    // +charged particle conditions
    for(int k=0; k<6; k++) {
            int i = ceil(ranfloat(nelem));
            int j = ceil(ranfloat(nelem));
            rho[i + j*nelem] = p[i + j*nelem] = 1;
    }

    // finite difference
    dx   = LENGTHX * 1 / (double) nelem;
    dy   = LENGTHY * 1 / (double) melem;
    dx2  = pow(dx, 2);
    dy2  = pow(dy, 2);

}

void Poisson2D::smooth(int l, int nu, double *v, double *f) {
    /* 
     * jacobi relaxation smoothy algorithm. iteratively solve \nabla^2 u = f
     * assumptions: (1) v, f are square lattices with 2^l+1 elements on each dimension
     *              (2) l is level number (index starts at +1)
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
    int n=pow(2,l); //nelem;
    int m=pow(2,l); //melem;
    int step=pow(2,nu);

    double trace = 0;
    double traceold = 0;
    
    for(int j=step; j<(m-1); j+=step) {
        for(int i=step; i<(n-1); i+=step) {
            // laplacian finite difference
            up          =v[i + (j+1)*n];
            down        =v[i + (j-1)*n];
            right       =v[(i+1) + j*n];
            left        =v[(i-1) + j*n];
            // u^n
            u[i + j*n]  =v[i + j*n];
            // u^n+1
            v[i + j*n]  =0.5*( dx2*(up + down)/(dx2 + dy2) + 
                                 dy2*(right + left)/(dx2 + dy2) +
                                 dx2*dy2*f[i + j*n]/(dx2 + dy2) );
            // error
            if(i==j) {
                trace    +=v[i + j*n];      // new trace
                traceold +=u[i + j*n];      // old trace
            }

        }
    }
    error = (pow(trace,2) - pow(traceold,2)) / pow(traceold,2);
   
    double val; 
    trace = 0;
    traceold = 0;
    for(int j=step; j<(m-1); j+=step) {
        for(int i=step; i<(n-1); i+=step) {
            // laplacian finite difference
            up          =utmp[i + (j+1)*n];
            down        =utmp[i + (j-1)*n];
            right       =utmp[(i+1) + j*n];
            left        =utmp[(i-1) + j*n];
            // u^n
            val         =utmp[i + j*n];
            // u^n+1
            utmp[i + j*n]=0.5*( dx2*(up + down)/(dx2 + dy2) + 
                                dy2*(right + left)/(dx2 + dy2) +
                                dx2*dy2*f[i + j*n]/(dx2 + dy2) );
            // error
            if(i==j) {
                trace    +=utmp[i + j*n];   // new trace
                traceold +=val;             // old trace
            }

        }
    }
    errortmp = (pow(trace,2) - pow(traceold,2)) / pow(traceold,2);
    
    if(DBG) printf("error = %lf\n", error);
}

void Poisson2D::defect(int l, int nu, double *v, double *f, double* dl) {
    /*
     * defect calculation.
     *
     *      dl = Ll ul - fl
     *
     * determines minus defect
     */
    if (DBG) printf("\nPoisson2D::defect()\n");
    
    double up, down, right, left, onsite;
    int n=pow(2,l); //n-elements;
    int m=pow(2,l); //m-elements;
    int step=pow(2,nu);

    double d2i=step * LENGTHX * 1 / (double) n;
    
    for(int j=step; j<(m-1); j+=step) { // interior defect points
        for(int i=step; i<(n-1); i+=step) {
            up          =v[i + (j+1)*n]; // Ll ul - fl operation
            down        =v[i + (j-1)*n];
            right       =v[(i+1) + j*n];
            left        =v[(i-1) + j*n];
            onsite      =v[i + j*n];
            /*dl[i + j*n] =0.5 * ( dx2*(up + down)/(dx2 + dy2) +
                                 dy2*(right + left)/(dx2 + dy2) +
                                 (dx2*dy2)*f[i + j*n]/(dx2 + dy2) ) - f[i + j*n];*/
            
            // minus defect (since \nabla^2 u = -p)
            // this is actually the pos defect
            dl[i+j*n]   =-d2i*(right + left + up + down - 4*onsite) + f[i+j*n];
            // now its neg defect
            dl[i+j*n]   =-dl[i+j*n];
        }
    }

    for(int i=0; i<n; i++) dl[i + 0]   = dl[i + (m-1)*n] = 0; // boundary defect points
    for(int j=0; j<m; j++) dl[0 + j*n] = dl[(n-1) + j*n] = 0;

}

void Poisson2D::restriction(int l, int nu, double *dl, double *rdl) {
    /*
     * restriction calculation.
     *
     *      dl-1 = r dl
     *
     * determines solution on l-1 level.
     */
    if (DBG) printf("\nPoisson2D::restriction()\n");
    
    double up, down, right, left, onsite;
    int n=pow(2,l); //n-elements;
    int m=pow(2,l); //m-elements;
    int step=pow(2,nu);

    for(int j=step; j<(m-1); j+=step) { // interior restrict points on l-1 level
        for(int i=step; i<(n-1); i+=step) {
            up          =dl[i + (j+1)*n]; // Ll ul - fl operation
            down        =dl[i + (j-1)*n];
            right       =dl[(i+1) + j*n];
            left        =dl[(i-1) + j*n];
            onsite      =dl[i + j*n];
            rdl[i + j*n] =0.125*(right + left + up + down) + 0.5*onsite;
        }
    }

}

void Poisson2D::prolongation(int l, int nu, double *rdl, double *dl, double* vtmp) {
    /*
     * restriction calculation.
     *
     *       \bar{u}l = pLl-1^-1[r dl]
     *
     * interpolates solution for the l'th level
     */
    if (DBG) printf("\nPoisson2D::restriction()\n");
    int n=pow(2,l); //n-elements;
    int m=pow(2,l); //m-elements;
    int step=pow(2,nu+1);

    // interpolation scheme
    for(int j=step; j<(m-1); j+=step) { // interior points on l-1 level
        for(int i=step; i<(n-1); i+=step) {
            dl[i + j*n] = rdl[i + j*n];
        }
    }

    // interpolate even numbered columns
    // and vertically odd rows.
    for(int j=1; j<(m-1); j+=step) {
        for(int i=0; i<(n-1); i+=step) {
            dl[i + j*n] = 0.5*(rdl[i + (j+1)*n] + rdl[i + (j-1)*n]);
        }
    }

    // interpolate odd numbered columns
    // and each row, interpolating horizontally.
    for(int j=0; j<(m-1); j+=1) {
        for(int i=1; i<(n-1); i+=step) {
            dl[i + j*n] = 0.5*(rdl[(i+1) + j*n] + rdl[(i-1) + j*n]);
        }
    }

    // add correction
    // ul^j+1 = ul^j - pvl-1
    double trace = 0;
    double traceold = 0;
    double vh;
    for(int j=1; j<(m-1); j+=1) {
        for(int i=1; i<(n-1); i+=1) {
            // u correction
            utmp[i + j*n]=vtmp[i + j*n] - dl[i + j*n];
        }
    }

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
        outfile << left << "u(x,y)\n";
        //outfile << left << "error" << "\n";
        outfile << std::scientific;
        outfile.precision(4);
        for(int j=1; j<(m-1); j++) {
            for(int i=1; i<(n-1); i++) {
                outfile << left << setw(12) << i*dx;
                outfile << left << setw(12) << j*dy;
                outfile << left << u[i + j*n] << "\n";
                //outfile << left << error << "\n";
            }
            outfile << "\n";
        }
    }
    else if (mode=="potentialtmp") {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# x";
        outfile << left << setw(12) << "y";
        outfile << left << "u(x,y)\n";
        //outfile << left << "error" << "\n";
        outfile << std::scientific;
        outfile.precision(4);
        for(int j=0; j<m; j++) {
            for(int i=0; i<n; i++) {
                outfile << left << setw(12) << i*dx;
                outfile << left << setw(12) << j*dy;
                outfile << left << utmp[i + j*n] << "\n";
                //outfile << left << errortmp << "\n";
            }
            outfile << "\n";
        }
    }
    else if (mode=="error") {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# index";
        outfile << left << "error\n";
        outfile << std::scientific;
        outfile.precision(4);
        outfile << left << setw(12) << index;
        outfile << left << error << "\n";
    }
    else if (mode=="errortmp") {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# index";
        outfile << left << "error_tmp\n";
        outfile << std::scientific;
        outfile.precision(4);
        outfile << left << setw(12) << index;
        outfile << left << errortmp << "\n";
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
                outfile << left << setw(12) << i*dx;
                outfile << left << setw(12) << j*dy;
                outfile << left << p[i + j*n] << "\n";
            }
            outfile << "\n";
        }
    }

}

// EXTRA STUFF


void Poisson2D::copy(int l, double* u, double* utmp) {
    int n=pow(2,l); //nelem;
    int m=pow(2,l); //melem;

    // internal points
    for(int j=1; j<(m-1); j++) {
        for(int i=1; i<(n-1); i++) {
            utmp[i + j*n] = u[i + j*n];
        }
    }

    // boundary conditions
    for(int i=0; i<n; i++) {
        utmp[i + 0*n] = 0;      // bottom bound
        utmp[0 + i*n] = 0;      // left bound
        utmp[i + (n-1)*n] = 0;  // top bound
        utmp[(n-1) + i*n] = 0;  // right bound
    }

}
