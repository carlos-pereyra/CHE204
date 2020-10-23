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

float* ComputeMorsePotential(long natoms, float *dr, float *ep);
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
    float *dxvec = (float *) malloc(sizeof(float)*natoms);       // size (natom by natom)
    float *dyvec = (float *) malloc(sizeof(float)*natoms);       // size (natom by natom)
    float *dzvec = (float *) malloc(sizeof(float)*natoms);       // size (natom by natom)
    float *drvec = (float *) malloc(sizeof(float)*natoms);       // size (natom by natom)
    //& jacobian force matrix 
    float *fxvec_old = (float *) malloc(sizeof(float)*natoms);   // size (natom by natom)
    float *fyvec_old = (float *) malloc(sizeof(float)*natoms);   // size (natom by natom)
    float *fzvec_old = (float *) malloc(sizeof(float)*natoms);   // size (natom by natom)
    float *fxvec = (float *) malloc(sizeof(float)*natoms);   // size (natom by natom)
    float *fyvec = (float *) malloc(sizeof(float)*natoms);   // size (natom by natom)
    float *fzvec = (float *) malloc(sizeof(float)*natoms);   // size (natom by natom)
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
        printf("\n=======================\n");
        // save current values into array (for old value reference). 
        MatrixCopy(natoms, ndim,   xyzmat, xyzmat_old);
        MatrixCopy(natoms, 1, fxvec,  fxvec_old);
        MatrixCopy(natoms, 1, fyvec,  fyvec_old);
        MatrixCopy(natoms, 1, fzvec,  fzvec_old);

        //printf("\n__XYZ MATRIX__\n");
        //PrintMatrix(natoms, natoms, xyzmat);
        
        // clear all forces (n x n)
        ClearMatrix(natoms, 1, fxvec);
        ClearMatrix(natoms, 1, fyvec);
        ClearMatrix(natoms, 1, fzvec);

        // compute dx, dy, dz, dr (n x n)
        ComputeDX(natoms, xyzmat, dxvec);
        ComputeDY(natoms, xyzmat, dyvec);
        ComputeDZ(natoms, xyzmat, dzvec);
        ComputeDR(natoms, xyzmat, drvec);
   
        printf("\n__DX_VEC__\n");
        PrintMatrix(natoms, 1, dxvec);
        printf("\n__DY_VEC__\n");
        PrintMatrix(natoms, 1, dyvec);
        printf("\n__DZ_VEC__\n");
        PrintMatrix(natoms, 1, dzvec);
        printf("\n__DR_VEC__\n");
        PrintMatrix(natoms, 1, drvec);

        // compute fx, fy, fz jacobian matrices (diagonal elements are atom forces) (n x n)
        ComputeMorseForceX(natoms, dxvec, drvec, fxvec);
        ComputeMorseForceY(natoms, dyvec, drvec, fyvec);
        ComputeMorseForceZ(natoms, dzvec, drvec, fzvec);

        printf("\n__FX_VECTOR__\n");
        PrintMatrix(natoms, 1, fxvec);
        printf("\n__FY_VECTOR__\n");
        PrintMatrix(natoms, 1, fyvec);
        printf("\n__FZ_VECTOR__\n");
        PrintMatrix(natoms, 1, fzvec);
        
        // compute new positions (returns new xyzmat)
        ComputeNewX(natoms, xyzmat, gammax, fxvec); 
        ComputeNewY(natoms, xyzmat, gammay, fyvec);
        ComputeNewZ(natoms, xyzmat, gammaz, fzvec);

        // compute optimized gamma ceof.
        gammax = ComputeNewGammaX(natoms, xyzmat, xyzmat_old, gammax, fxvec, fxvec_old);
        gammay = ComputeNewGammaY(natoms, xyzmat, xyzmat_old, gammay, fyvec, fyvec_old);
        gammaz = ComputeNewGammaZ(natoms, xyzmat, xyzmat_old, gammaz, fzvec, fzvec_old);

        // compute maximum fx, fy, fz
        int idX = GetAbsMaxElementIndex(natoms, 1, fxvec, 1); // stride is 1
        int idY = GetAbsMaxElementIndex(natoms, 1, fyvec, 1);
        int idZ = GetAbsMaxElementIndex(natoms, 1, fzvec, 1);

        // show absolute maximum element
        if(abs(fxvec[idX]) > abs(fyvec[idY]) && 
            abs(fxvec[idX]) > abs(fzvec[idZ])) {
            cout << std::scientific;
            cout.precision(2);
            cout << "\nMax force component is fx[" << idX;
            cout << "] = " << fxvec[idX] << "\n";
        }
        if(abs(fyvec[idY]) > abs(fxvec[idX]) && 
            abs(fyvec[idY]) > abs(fzvec[idZ])) {
            cout << std::scientific;
            cout.precision(2);
            cout << "\nMax force component is fy[" << idY;
            cout << "] = " << fyvec[idY] << "\n";
        }
        if(abs(fzvec[idZ]) > abs(fxvec[idX]) && 
            abs(fzvec[idZ]) > abs(fzvec[idY])) {
            cout << std::scientific;
            cout.precision(2);
            cout << "\nMax force component is fz[" << idZ;
            cout << "] = " << fyvec[idZ] << "\n";
        }
        if(abs(fxvec[idX])<(1e-4) || abs(fyvec[idY])<(1e-4) || abs(fzvec[idZ])<(1e-4) ) {
            break;
        }
        Coords2XYZFile(natoms, xyzmat, n);
    }

    // free memory 
    free(xyzmat); free(xyzmat_old);
    free(dxvec);     free(dyvec);     free(dzvec); free(drvec); 
    free(fxvec);     free(fyvec);     free(fzvec); 
    free(fxvec_old); free(fyvec_old); free(fzvec_old); 
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
        if( !((i+1)%m) ) { cout << "\n"; }
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

float* ComputeDX(long natoms, float *xyzmat, float *dxvec) {
    // provide natoms for the shape of the following matrices:
    //      - access (natoms x 3) xyzmat matrix
    //      - access (natoms x natoms) xyzmat
    // return dxmat
    if(DBG) printf("\nComputeDX()\n");
    for(long i = 0; i < natoms; i++) {
        dxvec[i] = xyzmat[0 + 0] - xyzmat[3*i + 0];
    }
    return dxvec;
}

float* ComputeDY(long natoms, float *xyzmat, float *dyvec) {
    // provide natoms for the shape of the following matrices:
    //      - access (natoms x 3) xyzmat matrix
    //      - access (natoms x natoms) dxmat
    // return dymat
    if(DBG) printf("\nComputeDY()\n");
    for(long i = 0; i < natoms; i++) {
        dyvec[i] = xyzmat[0 + 1] - xyzmat[3*i + 1];
    }
    return dyvec;
}

float* ComputeDZ(long natoms, float *xyzmat, float *dzvec) {
    // provide natoms for the shape of the following matrices:
    //      - access (natoms x 3) xyzmat matrix
    //      - access (natoms x 1) dzvec vector
    // return dzmat
    if(DBG) printf("\nComputeDZ()\n");
    for(long i = 0; i < natoms; i++) {
        dzvec[i] = xyzmat[0 + 2] - xyzmat[3*i + 2];
    }
    return dzvec;
}

float* ComputeDR(long natoms, float *xyzmat, float *drvec) {
    // provide natoms for the shape of the following matrices:
    //      - expects (natoms x 3) xyzmat matrix
    //      - expects (natoms x 1) drvec vector
    // return drmat
    float dx, dy, dz;
    if(DBG) printf("\nComputeDR()\n");
    for(long i = 0; i<natoms; i++) {
        dx = xyzmat[0 + 0] - xyzmat[i + 0];
        dy = xyzmat[0 + 1] - xyzmat[i + 1];
        dz = xyzmat[0 + 2] - xyzmat[i + 2];
        drvec[i] = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
    }
    return drvec;
}

// POTENTIAL CALCULATION

float* ComputeMorsePotential(long natoms, float *dr, float *ep) {
    // this function computes the morse potential. 
    // this function expects vectors *dr, *ep with length natoms
    //
    // return a list containing potential energy.
    float d = 2; float a = 1; float ro = 1.2; float ep_value; float ep_sum;
    long i, j;
    for(i = 0; i < natoms; i++) {
        ep_sum = 0;
        for(j = i + 1; j < natoms; j++) {
            ep_value = d * pow((1 - exp(-a*(dr[i]-ro))), 2);
            ep[j] += ep_value;
            ep_sum += ep_value;
        }
        ep[i] = ep_sum;
    }
    return ep;
}

// FORCE CALCULATIONS

float* ComputeMorseForceX(long natoms, float *dx, float *dr, float *fx) {
    // provide natoms for the shape of the following matrices:
    //      - dx matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fx matrix size (natoms x natoms)
    // return fx
    if(DBG) printf("\nComputeMorseForceX()\n");
    float d = 2; float a = 1; float ro = 1.2; float val = 0; float force_sum;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        force_sum = 0;
        for(j = i + 1; j < natoms; j++) {
            val = 2*d*a * ( exp(-2*a*(dr[j]-ro)) - exp(-a*(dr[j]-ro)) ) * dx[j] / dr[j];
            printf("val = %f\n", val);
            fx[j] -= val;
            force_sum += val;
            //fx[i*natoms + j] = val; 
            //fx[j*(natoms + 1)] -= val;
        }
        // particle force elements
        fx[i] += force_sum;
    }
    // diagonal forces in matrix
    /*for(long i = 0; i < natoms; i++) {
        for(long j = i+1; j < natoms; j++) {
            fx[i*(natoms + 1)] += fx[i*natoms + j];
        }
    }*/
    return fx;
}

float* ComputeMorseForceY(long natoms, float *dy, float *dr, float *fy) {
    // provide natoms for the shape of the following matrices:
    //      - dy matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fx matrix size (natoms x natoms)
    // return fy matrix
    if(DBG) printf("\nComputeMorseForceY()\n");
    float d = 2; float a = 1; float ro = 1.2; float val = 0; float force_sum;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        for(j = i + 1; j < natoms; j++) {
            val = 2*d*a * ( exp(-2*a*(dr[j]-ro)) - exp(-a*(dr[j]-ro)) ) * dy[j] / dr[j];
            fy[j] -= val;
            force_sum += val;
            //fy[i*natoms + j] = val;
            //fy[j*(natoms + 1)] -= val;
        }
        fy[i] += force_sum;
    }
    // diagonal forces in matrix
    /*for(long i = 0; i < natoms; i++) {
        for(long j = i+1; j < natoms; j++) {
            fy[i*(natoms + 1)] += fy[i*natoms + j];
        }
    }*/
    return fy;
}

float* ComputeMorseForceZ(long natoms, float *dz, float *dr, float *fz) {
    // provide natoms for the shape of the following matrices:
    //      - dz matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fz matrix size (natoms x natoms)
    // return fz matrix
    if(DBG) printf("\nComputeMorseForceZ()\n");
    float d = 2; float a = 1; float ro = 1.2; float val = 0; float force_sum;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        for(j = i + 1; j < natoms; j++) {
            val = 2*d*a * ( exp(-2*a*(dr[j]-ro)) - exp(-a*(dr[j]-ro)) ) * dz[j] / dr[j];
            fz[j] -= val;
            force_sum += val;
            //fz[i*natoms + j] = val;
            //fz[j*(natoms + 1)] -= val;
        }
        fz[i] += force_sum;
    }
    // diagonal forces in matrix
    /*for(long i = 0; i < natoms; i++) {
        for(long j = i+1; j < natoms; j++) {
            fz[i*(natoms + 1)] += fz[i*natoms + j];
        }
    }*/
    return fz;
}

// NEW POSITIONS

float* ComputeNewX(long n, float* xyzmat, float gx, float* fxvec) {
    // compute new x coorindate with new gamma and force computations
    //      eqn: compute x^(k+1)[i] = x[i] + gammax[i]*fx[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewX()\n");
    
    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 0] = xyzmat[i*ndim + 0] + gx * fxvec[i];
    }
    return xyzmat;
}

float* ComputeNewY(long n, float* xyzmat, float gy, float* fyvec) {
    // compute new y coorindate with new gamma and force computations
    //      eqn: compute y^(k+1)[i] = y[i] + gy * fy[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewY()\n");
    
    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 1] = xyzmat[i*ndim + 1] + gy * fyvec[i];
    }
    return xyzmat;
}

float* ComputeNewZ(long n, float* xyzmat, float gz, float* fzvec) {
    // compute new z coorindate with new gamma and force computations
    //      eqn: compute z^(k+1)[i] = z[i] + gz * fz[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewZ()\n");

    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 2] = xyzmat[i*ndim + 2] + gz * fzvec[i];
    }
    return xyzmat;
}

// NEW GAMMA

float ComputeNewGammaX(long n, float* xyz, float* xyzold, float gx, float* fx, float* fxold) {
    // compute new gamma for x direction
    //  expects xyz (natoms x 3) matrix
    //  expects xyzold (natoms x 3) matrix
    //  expects gx scalar
    //  expects fx (natoms x 1) vector
    //  expects fxold (natoms x 1) vector
    if(DBG) printf("\nComputeNewGammaX()\n");
    long i;
    float dx, df;
    float numerator_sum = 0;
    float denominator_sum = 0;
    for(i = 0; i < n; i++) {
        dx = xyz[ i*3 + 0 ] - xyzold[ i*3 + 0 ];
        df = fx[i]   - fxold[i];
        numerator_sum += dx * df;
        denominator_sum += pow(df, 2);
        //printf("dx = %f df = %f\n", dx, df);
    }
    if((numerator_sum == 0) || (denominator_sum == 0)) return 0;
    gx = numerator_sum / denominator_sum;
    return gx;
}

float ComputeNewGammaY(long n, float* xyzmat, float* xyzmatold, float gy, float* fy, float* fyold) {
    // compute new gamma for y direction
    if(DBG) printf("\nComputeNewGammaY()\n");
    long i; 
    float dy, df;
    float numerator_sum = 0; 
    float denominator_sum = 0;
    for(i = 0; i < n; i++) { 
        dy = xyzmat[ i*3 + 1 ] - xyzmatold[ i*3 + 1 ];
        df = fy[i]   - fyold[i];
        numerator_sum += dy * df;
        denominator_sum += pow(df, 2);
        //printf("dy = %f df = %f\n", dy, df);
    }
    if((numerator_sum == 0) || (denominator_sum == 0)) return 0;
    gy = numerator_sum / denominator_sum;
    return gy;
}

float ComputeNewGammaZ(long n, float* xyzmat, float* xyzmatold, float gz, float* fz, float* fzold) {
    // compute new gamma for z direction
    if(DBG) printf("\nComputeNewGammaZ()\n");
    long i; 
    float dz, df;
    float numerator_sum = 0;
    float denominator_sum = 0;
    for(i = 0; i < n; i++) {
        dz = xyzmat[ i*3 + 2 ] - xyzmatold[ i*3 + 2 ];
        df = fz[i]   - fzold[i];
        numerator_sum += dz * df;
        denominator_sum += pow(df, 2);
        //printf("dz = %f df = %f\n", dz, df);
    }
    if((numerator_sum == 0) || (denominator_sum == 0)) return 0;
    gz = numerator_sum / denominator_sum;
    return gz;
}
