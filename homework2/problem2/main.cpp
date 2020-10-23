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
#include "MersenneTwister.h"

#define DBG 0

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

// (on OS X) compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -framework Accelerate -o test main.cpp
//
//      ./test
//

#define ranfloat(w) ( w * ((rand() / (double) RAND_MAX)  - 0.5) )
using namespace std;

float* readInput(string filename, long nlines, float *mat);
float* ClearMatrix(long n, long m, float* mat);
float* SetMatrixElement(long n, long m, float *mat, long index, long value);
int GetAbsMaxElementIndex(long n, long m, float* mat, long skip);

// GDA RELATED
float* SetInitialGuessPositions(long n, long m, MTRand *, float* xyz, float* newxyz);
float* ComputeGammaX(long n, long m, float* xyzmat, float* xyzmat0, float* gamma);
float* MatrixCopy(long n, long m, float* x, float* y);
float* VectorExtraction(long n, long m, float*, long, float*);
// EXTRA
float* setSMatrix(long n, long m, float sxy, float syz, float szy); // ignore
float* setXLatticePosition(long n, float *x, float dx);             // ignore
float* setYLatticePosition(long n, float *y, float dy);             // ignore
float* setZLatticePosition(long n, float *z, float dz);             // ignore

// IO
void PrintMatrix(const long n, const long m, float *);
void PrintMatrixDiag(const long n, const long m, float *);
void Coords2XYZFile(long n, float* xyz, long index);

// PHYSICS
float* MorsePotential(long s, long e, float *r, float *ep);

float* ComputeDX(long natoms, float *xyzmat, float *dxmat);
float* ComputeDY(long natoms, float *xyzmat, float *dymat);
float* ComputeDZ(long natoms, float *xyzmat, float *dzmat);
float* ComputeDR(long natoms, float *xyzmat, float *drmat);

float* ComputeMorseForceX(long n, float *dx, float *dr, float *fx);
float* ComputeMorseForceY(long n, float *dy, float *dr, float *fy);
float* ComputeMorseForceZ(long n, float *dz, float *dr, float *fz);

float* ComputeNewX(long n, float* xyzmat, float gx, float* fxmat);
float* ComputeNewY(long n, float* xyzmat, float gy, float* fymat);
float* ComputeNewZ(long n, float* xyzmat, float gz, float* fzmat);

float ComputeNewGammaX(long n, float* xyzmat, float* xyzmatold, float gx, float* fxmat, float* fxmatold);
float ComputeNewGammaY(long n, float* xyzmat, float* xyzmatold, float gy, float* fymat, float* fymatold);
float ComputeNewGammaZ(long n, float* xyzmat, float* xyzmatold, float gz, float* fzmat, float* fzmatold);

/*
float* ComputeGamma(long n, long m, float* xyzmat, float* xyzmat0, float* gamma);
float* ComputeGammaVec(long n, float* r, float* r0, float* fr, float* fr0, float* gamma);
*/

int main(int argc, char** argv) {
    // setup and initialization
    long natoms = 3; long ndim = 3; // make sure natoms matches the number of atoms in .xyz file
    long numruns = 3;
    // MATRIX
    //
    //& coords
    float *xyzmat_old = (float *) malloc(sizeof(float)*natoms*ndim);    // size (natom by ndim)
    float *xyzmat = (float *) malloc(sizeof(float)*natoms*ndim);    // size (natom by ndim)
    //& pair distance
    float *dxmat = (float *) malloc(sizeof(float)*natoms*natoms);       // size (natom by natom)
    float *dymat = (float *) malloc(sizeof(float)*natoms*natoms);       // size (natom by natom)
    float *dzmat = (float *) malloc(sizeof(float)*natoms*natoms);       // size (natom by natom)
    float *drmat = (float *) malloc(sizeof(float)*natoms*natoms);       // size (natom by natom)
    //& jacobian force matrix 
    float *fxmat_old = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *fymat_old = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *fzmat_old = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *fxmat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *fymat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *fzmat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    // maximum force array
    float *maxf = (float *) malloc(sizeof(float)*3);   // size (natom by natom)
    // random number generator
    unsigned long int now = static_cast<unsigned long int>(time(NULL));
    MTRand *mtrand = new MTRand(now);

    // read xyz positions and 
    readInput("coordinates/Test.xyz", natoms, xyzmat);

    // initial guess gamma values 
    float gammax = 1; //mtrand->randNorm(0,1);
    float gammay = 1; //mtrand->randNorm(0,1);
    float gammaz = 1; //mtrand->randNorm(0,1);
    
    for(long n=0; n<numruns; n++) {
        printf("\n=======================");
        printf("\n__ iteration = %ld __\n", n);   
        printf("\n=======================");
        // save current values into array (for old value reference). 
        MatrixCopy(natoms, ndim,   xyzmat, xyzmat_old);
        MatrixCopy(natoms, natoms, fxmat,  fxmat_old);
        MatrixCopy(natoms, natoms, fymat,  fymat_old);
        MatrixCopy(natoms, natoms, fzmat,  fzmat_old);

        //printf("\n__XYZ MATRIX__\n");
        //PrintMatrix(natoms, natoms, xyzmat);
        
        // clear all forces (n x n)
        ClearMatrix(natoms, natoms, fxmat);
        ClearMatrix(natoms, natoms, fymat);
        ClearMatrix(natoms, natoms, fzmat);

        /*
        printf("\n__FX MATRIX__\n");
        PrintMatrix(natoms, natoms, fxmat);
        printf("\n__FY MATRIX__\n");
        PrintMatrix(natoms, natoms, fymat);
        printf("\n__FZ MATRIX__\n");
        PrintMatrix(natoms, natoms, fzmat);
        */

        // compute dx, dy, dz, dr (n x n)
        ComputeDX(natoms, xyzmat, dxmat);
        ComputeDY(natoms, xyzmat, dymat);
        ComputeDZ(natoms, xyzmat, dzmat);
        ComputeDR(natoms, xyzmat, drmat);
   
        printf("\n__DX_MATRIX__\n");
        PrintMatrix(natoms, natoms, dxmat);
        printf("\n__DY_MATRIX__\n");
        PrintMatrix(natoms, natoms, dymat);
        printf("\n__DZ_MATRIX__\n");
        PrintMatrix(natoms, natoms, dzmat);

        // compute fx, fy, fz jacobian matrices (diagonal elements are atom forces) (n x n)
        ComputeMorseForceX(natoms, dxmat, drmat, fxmat);
        ComputeMorseForceY(natoms, dymat, drmat, fymat);
        ComputeMorseForceZ(natoms, dzmat, drmat, fzmat);

        // compute new positions (returns new xyzmat)
        ComputeNewX(natoms, xyzmat, gammax, fxmat); 
        ComputeNewY(natoms, xyzmat, gammay, fymat);
        ComputeNewZ(natoms, xyzmat, gammaz, fzmat);

        printf("\n__FX'_MATRIX__\n");
        PrintMatrix(natoms, natoms, fxmat);
        printf("\n__FY'_MATRIX__\n");
        PrintMatrix(natoms, natoms, fymat);
        printf("\n__FZ'_MATRIX__\n");
        PrintMatrix(natoms, natoms, fzmat);

        // compute optimized gamma ceof.
        gammax = ComputeNewGammaX(natoms, xyzmat, xyzmat_old, gammax, fxmat, fxmat_old);
        gammay = ComputeNewGammaY(natoms, xyzmat, xyzmat_old, gammay, fymat, fymat_old);
        gammaz = ComputeNewGammaZ(natoms, xyzmat, xyzmat_old, gammaz, fzmat, fzmat_old);

        // compute maximum fx, fy, fz
        int max_idx_x = GetAbsMaxElementIndex(natoms, natoms, fxmat, 1); // stride is 1
        int max_idx_y = GetAbsMaxElementIndex(natoms, natoms, fymat, 1);
        int max_idx_z = GetAbsMaxElementIndex(natoms, natoms, fzmat, 1);

        // show absolute maximum element
        if(abs(fxmat[max_idx_x]) > abs(fymat[max_idx_y]) && 
            abs(fxmat[max_idx_x]) > abs(fzmat[max_idx_z])) {
            cout << std::scientific;
            cout.precision(2);
            cout << "\nMax force component is fx[" << max_idx_x;
            cout << "] = " << fxmat[max_idx_x] << "\n";
        }
        if(abs(fymat[max_idx_y]) > abs(fxmat[max_idx_x]) && 
            abs(fymat[max_idx_y]) > abs(fzmat[max_idx_z])) {
            cout << std::scientific;
            cout.precision(2);
            cout << "\nMax force component is fy[" << max_idx_y;
            cout << "] = " << fymat[max_idx_y] << "\n";
        }
        if(abs(fzmat[max_idx_z]) > abs(fxmat[max_idx_x]) && 
            abs(fzmat[max_idx_z]) > abs(fzmat[max_idx_y])) {
            cout << std::scientific;
            cout.precision(2);
            cout << "\nMax force component is fz[" << max_idx_z;
            cout << "] = " << fymat[max_idx_z] << "\n";
        }
        if(abs(fxmat[max_idx_x])<(1e-4) || abs(fymat[max_idx_y])<(1e-4) || abs(fzmat[max_idx_z])<(1e-4) ) {
            break;
        }
        Coords2XYZFile(natoms, xyzmat, n);
    }

    // free memory 
    free(xyzmat); free(xyzmat_old);
    free(dxmat);     free(dymat);     free(dzmat); free(drmat); 
    free(fxmat);     free(fymat);     free(fzmat); 
    free(fxmat_old); free(fymat_old); free(fzmat_old); 
    free(maxf);
    delete mtrand;
}

// ========// ========// ========// ========// =============

float* readInput(string filename, long nlines, float *mat) {
    // read in positions from .xyz file
    ifstream infile(filename);
    string elem, line;
    float x, y, z;
    long counter, natoms;
    
    if(!infile.is_open()) exit(1);
    
    counter=0;
    std::getline(infile, line); std::istringstream iss(line); iss >> natoms; // line 1
    std::cout << "natoms = " << natoms << "\n"; 
    std::getline(infile, line);  // line 2
    while(std::getline(infile, line) || counter<nlines) {
        // read x y z coords
        std::istringstream iss(line);

        if (iss >> elem >> x >> y >> z) cout << counter << " elem = " << elem << " x = " << x << " y = " << y << " z " << z << "\n";

        mat[3*counter] = x;
        mat[3*counter+1] = y;
        mat[3*counter+2] = z;

        counter++;
    }
    infile.close();
    return mat;
}

float* ClearMatrix(long n, long m, float *vec) {
    // clear vector that is a real n by m matrix
    if(DBG) printf("\nClearMatrix()\n");
    if(__APPLE__) {
        // catlas_sset(int <numelem>, int <val>, float* <vec>, incx <stride>)
        catlas_sset(n*m, 0, vec, 1);
    } else {
        for(long i = 0; i < n; i++) {
            for(long j = 0; j < m; j++) {
                vec[i*n + j] = 0;
            }
        }
    }
    return vec;
}

float* SetMatrixElement(long n, long m, float *mat, long index, long value) {
    // set vector element (matrix is a real n by m) to some value.
    mat[index] = value;
    return mat;
}

int GetAbsMaxElementIndex(long n, long m, float *vec, long skip) {
    // get abs(max) element in vector that is a real n by m matrix
    // skip is the stride constant
    int index; float max;
    if(__APPLE__) {
        index = cblas_isamax(n*m, vec, skip);
    } else {
        index = 0; max = 0;
        for(long i = 0; i < n; i++) {
            for(long j = 0; j < m; j++) {
                if( max < abs(vec[i*n + j]) ) max = abs(vec[i*n + j]);
            }
        }
    }
    return index;
}

// GDA RELATED
float* SetInitialGuessPositions(long n, long m, MTRand* mtrand, float* xyz, float* newxyz) {
    // setup initial guess positions
    // returns new xyz real matrix that is size n by m
    for(long i = 0; i < n; i++) {
        for(long j = 0; j < m; j++) {
            newxyz[i*m + j] = xyz[i*m + j] + mtrand->randNorm(0,1);
        }
    }
    return newxyz;
}

float* MatrixCopy(long n, long m, float* x, float* y) {
    // copy x matrix to y matrix
    if(DBG) printf("\nMatrixCopy()\n");
    if(__APPLE__) {
        cblas_scopy(n*m, x, 1, y, 1);
    }
    return y;
}

float* VectorExtraction(long n, long m, float* mat, long stride, float* vec) {
    // extract column vector from matrix
    if(__APPLE__) {
        cblas_scopy(n, mat, stride, vec, 1);
    }
    return vec;
}

// IO
void PrintMatrix(const long n, const long m, float *vec) {
    // print matrix
    for(long i = 0; i < n*m; i++) {
        //cout << std::scientific;
        cout.precision(2);
        cout << setw(8) << vec[i] << "\t";
        if( !((i+1)%n) ) { cout << "\n"; }
    }
}

void PrintMatrixDiag(const long n, const long m, float *vec) {
    // print matrix diagonal (assumes vec is real and square)
    for(long i = 0, j = 0; j < n*m; i++, j = i*(n+1)) {
        cout << std::scientific;
        cout.precision(2);
        cout << setw(8) << vec[j] << setw(8) << "i = " << j << "\n";
    }
}

void Coords2XYZFile(long n, float *xyz, long index) {
    // write coordinates to xyz file
    string filename_xyz = "data/output.xyz";
    ofstream myfile;
    long i;
    // open new file 
    if (index == 0) myfile.open(filename_xyz, std::fstream::out);
    // append to file
    else myfile.open(filename_xyz, std::fstream::app);

    myfile << n << "\n";
    myfile << "(n=" << index << ")!\n";
    for (i=0; i<n; i++) {
        // write atom info to file.
        myfile << "Na" << "\t";
        myfile << xyz[i*3+0] << "\t";
        myfile << xyz[i*3+1] << "\t";
        myfile << xyz[i*3+2] << "\n";
    }
    myfile.close();
}

// PHYSICS

float* MorsePotential(long s, long e, float *r, float *ep) {
    // this function computes the morse potential.
    //
    // return a list containing potential energy.
    float d = 2;
    float a = 1;
    float ro = 1.2;
    for(long i=s; i<e; i++){
        ep[i] = d*pow((1 - exp(-a*(r[i]-ro))), 2);
            printf("%ld -> %ld dr = %f ep = %f\n", s, i, r[i], ep[i]);
    }
    return ep;
}

// SEPARATION DISTANCES

float* ComputeDX(long natoms, float *xyzmat, float *dxmat) {
    // provide natoms for the shape of the following matrices:
    //      - access (natoms x 3) xyzmat matrix
    //      - access (natoms x natoms) xyzmat
    // return dxmat
    if(DBG) printf("\nComputeDX()\n");
    for(long i=0; i<natoms; i++) {
        for(long j=i+1; j<natoms; j++) {
            dxmat[natoms*i + j] = xyzmat[3*i] - xyzmat[3*j];
        }
    }
    return dxmat;
}

float* ComputeDY(long natoms, float *xyzmat, float *dymat) {
    // provide natoms for the shape of the following matrices:
    //      - access (natoms x 3) xyzmat matrix
    //      - access (natoms x natoms) dxmat
    // return dymat
    if(DBG) printf("\nComputeDY()\n");
    for(long i=0; i<natoms; i++) {
        for(long j=i+1; j<natoms; j++) {
            dymat[natoms*i + j] = xyzmat[3*i + 1] - xyzmat[3*j + 1];
        }
    }
    return dymat;
}

float* ComputeDZ(long natoms, float *xyzmat, float *dzmat) {
    // provide natoms for the shape of the following matrices:
    //      - access (natoms x 3) xyzmat matrix
    //      - access (natoms x natoms) dxmat
    // return dzmat
    if(DBG) printf("\nComputeDZ()\n");
    for(long i=0; i<natoms; i++) {
        for(long j=i+1; j<natoms; j++) {
            dzmat[natoms*i + j] = xyzmat[3*i + 2] - xyzmat[3*j + 2];
        }
    }
    return dzmat;
}

float* ComputeDR(long natoms, float *xyzmat, float *drmat) {
    // provide natoms for the shape of the following matrices:
    //      - access (natoms x 3) xyzmat matrix
    //      - access (natoms x natoms) dxmat
    // return drmat
    if(DBG) printf("\nComputeDR()\n");
    for(long i=0; i<natoms; i++) {
        for(long j=i+1; j<natoms; j++) {
            float dx = xyzmat[3*i] - xyzmat[3*j];
            float dy = xyzmat[3*i+1] - xyzmat[3*j+1];
            float dz = xyzmat[3*i+2] - xyzmat[3*j+2];
            drmat[natoms*i + j] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
        }
    }
    return drmat;
}

// FORCE CALCULATIONS

float* ComputeMorseForceX(long natoms, float *dx, float *dr, float *fx) {
    // provide natoms for the shape of the following matrices:
    //      - dx matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fx matrix size (natoms x natoms)
    // return fx
    if(DBG) printf("\nComputeMorseForceX()\n");
    float d = 2; float a = 1; float ro = 1.2; float val = 0;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        for(j = i + 1; j < natoms; j++) {
            val = 2*d*a * ( exp(-2*a*(dr[i*natoms + j]-ro)) - exp(-a*(dr[i*natoms + j]-ro)) ) * dx[i*natoms + j] / dr[i*natoms + j];
            fx[i*natoms + j] = val; 
            fx[j*(natoms + 1)] -= val;
        }
    }
    // diagonal forces in matrix
    for(long i = 0; i < natoms; i++) {
        for(long j = i+1; j < natoms; j++) {
            fx[i*(natoms + 1)] += fx[i*natoms + j];
        }
    }
    return fx;
}

float* ComputeMorseForceY(long natoms, float *dy, float *dr, float *fy) {
    // provide natoms for the shape of the following matrices:
    //      - dy matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fx matrix size (natoms x natoms)
    // return fy matrix
    if(DBG) printf("\nComputeMorseForceY()\n");
    float d = 2; float a = 1; float ro = 1.2; float val = 0;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        for(j = i + 1; j < natoms; j++) {
            val = 2*d*a * ( exp(-2*a*(dr[i*natoms + j]-ro)) - exp(-a*(dr[i*natoms + j]-ro)) ) * dy[i*natoms + j] / dr[i*natoms + j];
            fy[i*natoms + j] = val;
            fy[j*(natoms + 1)] -= val;
        }
    }
    // diagonal forces in matrix
    for(long i = 0; i < natoms; i++) {
        for(long j = i+1; j < natoms; j++) {
            fy[i*(natoms + 1)] += fy[i*natoms + j];
        }
    }
    return fy;
}

float* ComputeMorseForceZ(long natoms, float *dz, float *dr, float *fz) {
    // provide natoms for the shape of the following matrices:
    //      - dz matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fz matrix size (natoms x natoms)
    // return fz matrix
    if(DBG) printf("\nComputeMorseForceZ()\n");
    float d = 2; float a = 1; float ro = 1.2; float val = 0;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        for(j = i + 1; j < natoms; j++) {
            val = 2*d*a * ( exp(-2*a*(dr[i*natoms + j]-ro)) - exp(-a*(dr[i*natoms + j]-ro)) ) * dz[i*natoms + j] / dr[i*natoms + j];
            fz[i*natoms + j] = val;
            fz[j*(natoms + 1)] -= val;
        }
    }
    // diagonal forces in matrix
    for(long i = 0; i < natoms; i++) {
        for(long j = i+1; j < natoms; j++) {
            fz[i*(natoms + 1)] += fz[i*natoms + j];
        }
    }
    return fz;
}

// NEW POSITIONS

float* ComputeNewX(long n, float* xyzmat, float gx, float* fxmat) {
    // compute new x coorindate with new gamma and force computations
    //      eqn: compute x^(k+1)[i] = x[i] + gammax[i]*fx[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewX()\n");
    
    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 0] = xyzmat[i*ndim + 0] + gx * fxmat[i*(n+1)];
    }
    return xyzmat;
}

float* ComputeNewY(long n, float* xyzmat, float gy, float* fymat) {
    // compute new y coorindate with new gamma and force computations
    //      eqn: compute y^(k+1)[i] = y[i] + gy * fy[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewY()\n");
    
    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 1] = xyzmat[i*ndim + 1] + gy * fymat[i*(n+1)];
    }
    return xyzmat;
}

float* ComputeNewZ(long n, float* xyzmat, float gz, float* fzmat) {
    // compute new z coorindate with new gamma and force computations
    //      eqn: compute z^(k+1)[i] = z[i] + gz * fz[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewZ()\n");

    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 2] = xyzmat[i*ndim + 2] + gz * fzmat[i*(n+1)];
    }
    return xyzmat;
}

// NEW GAMMA

float ComputeNewGammaX(long n, float* xyzmat, float* xyzmatold, float gx, float* fxmat, float* fxmatold) {
    // compute new gamma for x direction
    if(DBG) printf("\nComputeNewGammaX()\n");
    long i;
    float dx, df;
    float numerator_sum = 0;
    float denominator_sum = 0;
    for(i = 0; i < n; i++) {
        dx = xyzmat[ i*3 + 0 ] - xyzmatold[ i*3 + 0 ];
        df = fxmat[ i*(n+1) ]   - fxmatold[ i*(n+1) ];
        numerator_sum += dx * df;
        denominator_sum += pow(df, 2);
        //printf("dx = %f df = %f\n", dx, df);
    }
    if((numerator_sum == 0) || (denominator_sum == 0)) return 0;
    gx = numerator_sum / denominator_sum;
    return gx;
}

float ComputeNewGammaY(long n, float* xyzmat, float* xyzmatold, float gy, float* fymat, float* fymatold) {
    // compute new gamma for y direction
    if(DBG) printf("\nComputeNewGammaY()\n");
    long i; 
    float dy, df;
    float numerator_sum = 0; 
    float denominator_sum = 0;
    for(i = 0; i < n; i++) { 
        dy = xyzmat[ i*3 + 1 ] - xyzmatold[ i*3 + 1 ];
        df = fymat[ i*(n+1) ]   - fymatold[ i*(n+1) ];
        numerator_sum += dy * df;
        denominator_sum += pow(df, 2);
        //printf("dy = %f df = %f\n", dy, df);
    }
    if((numerator_sum == 0) || (denominator_sum == 0)) return 0;
    gy = numerator_sum / denominator_sum;
    return gy;
}

float ComputeNewGammaZ(long n, float* xyzmat, float* xyzmatold, float gz, float* fzmat, float* fzmatold) {
    // compute new gamma for z direction
    if(DBG) printf("\nComputeNewGammaZ()\n");
    long i; 
    float dz, df;
    float numerator_sum = 0;
    float denominator_sum = 0;
    for(i = 0; i < n; i++) {
        dz = xyzmat[ i*3 + 2 ] - xyzmatold[ i*3 + 2 ];
        df = fzmat[ i*(n+1) ]   - fzmatold[ i*(n+1) ];
        numerator_sum += dz * df;
        denominator_sum += pow(df, 2);
        //printf("dz = %f df = %f\n", dz, df);
    }
    if((numerator_sum == 0) || (denominator_sum == 0)) return 0;
    gz = numerator_sum / denominator_sum;
    return gz;
}
