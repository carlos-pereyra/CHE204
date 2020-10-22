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
#define DBG 0
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
void Coords2XYZFile(const long n, float result, float mcresult);

// PHYSICS
float* MorsePotential(long s, long e, float *r, float *ep);
float* ComputeDX(long natoms, float *xyzmat, float *dxmat);
float* ComputeDY(long natoms, float *xyzmat, float *dymat);
float* ComputeDZ(long natoms, float *xyzmat, float *dzmat);
float* ComputeDR(long natoms, float *xyzmat, float *drmat);
float* ComputeMorseForceX(long natoms, float *dx, float *dr, float *fx);
float* ComputeMorseForceY(long natoms, float *dy, float *dr, float *fy);
float* ComputeMorseForceZ(long natoms, float *dz, float *dr, float *fz);

float* ComputeGamma(long n, long m, float* xyzmat, float* xyzmat0, float* gamma);
float* ComputeGammaVec(long n, float* r, float* r0, float* fr, float* fr0, float* gamma);

int main(int argc, char** argv) {
    // setup and initialization
    long natoms = 3; long ndim = 3; // make sure natoms matches the number of atoms in .xyz file
    long numruns = 1;
    // coords
    float *xyz0mat = (float *) malloc(sizeof(float)*natoms*ndim);   // size (natom by ndim)
    float *xyz1mat = (float *) malloc(sizeof(float)*natoms*ndim);   // size (natom by ndim)
    // MATRIX
    //
    //& pair distance
    float *dx0mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *dy0mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *dz0mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *dr0mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *dx1mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *dy1mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *dz1mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *dr1mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    //& jacobian force matrix 
    float *fx0mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *fy0mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *fz0mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *fx1mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *fy1mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    float *fz1mat = (float *) malloc(sizeof(float)*natoms*natoms);  // size (natom by natom)
    // VECTORS
    //
    //& position
    float *x0vec = (float *) malloc(sizeof(float)*natoms);          // size (natom by 1)
    float *y0vec = (float *) malloc(sizeof(float)*natoms);          // size (natom by 1)
    float *z0vec = (float *) malloc(sizeof(float)*natoms);          // size (natom by 1)
    float *x1vec = (float *) malloc(sizeof(float)*natoms);          // size (natom by 1)
    float *y1vec = (float *) malloc(sizeof(float)*natoms);          // size (natom by 1)
    float *z1vec = (float *) malloc(sizeof(float)*natoms);          // size (natom by 1)
    //& forces
    float *fx0vec = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    float *fy0vec = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    float *fz0vec = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    float *fx1vec = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    float *fy1vec = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    float *fz1vec = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    //& gamma
    float *gammax = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    float *gammay = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    float *gammaz = (float *) malloc(sizeof(float)*natoms);         // size (natom by 1)
    // maximum force array
    float *maxf = (float *) malloc(sizeof(float)*3);   // size (natom by natom)
    // random number generator
    unsigned long int now = static_cast<unsigned long int>(time(NULL));
    MTRand *mtrand = new MTRand(now);

    // read xyz positions and 
    readInput("coordinates/Test.xyz", natoms, xyz1mat);

    // initial guess values 
    SetInitialGuessPositions(natoms, ndim, mtrand, xyz1mat, xyz0mat);

    VectorExtraction(natoms, ndim, xyz0mat+0, ndim, x0vec);
    VectorExtraction(natoms, ndim, xyz0mat+1, ndim, y0vec);
    VectorExtraction(natoms, ndim, xyz0mat+2, ndim, z0vec);
   
    ClearMatrix(natoms, natoms, fx0mat);
    ClearMatrix(natoms, natoms, fy0mat);
    ClearMatrix(natoms, natoms, fz0mat);
    printf("\nFXYZ0__\n"); 
    PrintMatrix(natoms, natoms, fx0mat);

    // compute dx, dy, dz, dr (n x n)
    ComputeDX(natoms, xyz0mat, dx0mat);
    ComputeDY(natoms, xyz0mat, dy0mat);
    ComputeDZ(natoms, xyz0mat, dz0mat);
    ComputeDR(natoms, xyz0mat, dr0mat);

    // compute fx, fy, fz jacobian matrices (diagonal elements are atom forces) (n x n)
    ComputeMorseForceX(natoms, dx0mat, dr0mat, fx0mat);
    ComputeMorseForceY(natoms, dy0mat, dr0mat, fy0mat);
    ComputeMorseForceZ(natoms, dz0mat, dr0mat, fz0mat);
    
    // get fx, fy, fz components along diagonal
    VectorExtraction(natoms, natoms, fx0mat, natoms+1, fx0vec);
    VectorExtraction(natoms, natoms, fy0mat, natoms+1, fy0vec);
    VectorExtraction(natoms, natoms, fz0mat, natoms+1, fz0vec);
    /* 
    for(long n=0; n<numruns; n++) {    
        // clear all forces (n x n)
        ClearMatrix(natoms, natoms, fx1mat);
        ClearMatrix(natoms, natoms, fy1mat);
        ClearMatrix(natoms, natoms, fz1mat);
    
        // compute dx, dy, dz, dr (n x n)
        ComputeDX(natoms, xyz1mat, dx1mat);
        ComputeDY(natoms, xyz1mat, dy1mat);
        ComputeDZ(natoms, xyz1mat, dz1mat);
        ComputeDR(natoms, xyz1mat, dr1mat);
    
        // compute fx, fy, fz jacobian matrices (diagonal elements are atom forces) (n x n)
        ComputeMorseForceX(natoms, dxmat, drmat, fxmat);
        ComputeMorseForceY(natoms, dymat, drmat, fymat);
        ComputeMorseForceZ(natoms, dzmat, drmat, fzmat);
        
        // compute optimized gamma ceof.
        //ComputeGamma(natoms, ndim, xyzmat, xyzmat0, gamma); // matrix version
        ComputeGammaX(natoms, x1vec, x0vec, fx1vec, fx0vec, gammaxvec);
        ComputeGammaX(natoms, y1vec, y0vec, fy1vec, fy0vec, gammayvec);
        ComputeGammaX(natoms, z1vec, z0vec, fz1vec, fz0vec, gammazvec);

        // save current xyz coords to xyz-old
        VectorCopy(natoms, x1vec, x0vec);

        // compute x^(k+1)[i] = x[i] + gammax[i]*fx[i]
        // compute y^(k+1)[i] = y[i] + gammay[i]*fy[i]
        // compute z^(k+1)[i] = z[i] + gammaz[i]*fz[i]
        ComputeNewX(natoms, x1vec, gammaxvec, fx1vec); // returns x1vec
        ComputeNewY(natoms, y1vec, gammayvec, fy1vec); // returns y1vec
        ComputeNewZ(natoms, z1vec, gammazvec, fz1vec); // returns z1vec

        // compute maximum fx, fy, fz
        int max_idx_x = GetAbsMaxElementIndex(natoms, natoms, fxmat, 1); // stride is 1
        int max_idx_y = GetAbsMaxElementIndex(natoms, natoms, fymat, 1);
        int max_idx_z = GetAbsMaxElementIndex(natoms, natoms, fzmat, 1);
    
        // show absolute maximum element
        if(abs(fxmat[max_idx_x]) > abs(fymat[max_idx_y]) && abs(fxmat[max_idx_x]) > abs(fzmat[max_idx_z])) {
            printf("\nMax force component is fx[%d] = %f\n", max_idx_x, fxmat[max_idx_x]);
        }
        if(abs(fymat[max_idx_y]) > abs(fxmat[max_idx_x]) && abs(fymat[max_idx_y]) > abs(fzmat[max_idx_z])) {
            printf("\nMax force component is fy[%d] = %f\n", max_idx_y, fymat[max_idx_y]);
        }
        if(abs(fzmat[max_idx_z]) > abs(fxmat[max_idx_x]) && abs(fzmat[max_idx_z]) > abs(fzmat[max_idx_y])) {
            printf("\nMax force component is fz[%d] = %f\n", max_idx_z, fzmat[max_idx_z]);
        }
    }
    */
    // print summary
    printf("\nXYZ__\n\n");
    PrintMatrix(natoms, ndim, xyz1mat);
    printf("\nXYZ0__\n\n");
    PrintMatrix(natoms, ndim, xyz0mat);
    printf("\nX0_\n\n");
    PrintMatrix(natoms, 1, x0vec);
    printf("\nY0_\n\n");
    PrintMatrix(natoms, 1, y0vec);
    printf("\nZ0_\n\n");
    PrintMatrix(natoms, 1, z0vec);

    printf("\nFX0_\n\n");
    PrintMatrix(natoms, 1, fx0vec);
    // free XYZ
    free(xyz0mat); free(xyz1mat);
    // free MATRIX
    free(dx0mat); free(dy0mat); free(dz0mat); free(dr0mat); 
    free(dx1mat); free(dy1mat); free(dz1mat); free(dr1mat); 
    free(fx0mat); free(fy0mat); free(fz0mat); 
    free(fx1mat); free(fy1mat); free(fz1mat); 
    // free VECTORS
    free(x0vec);  free(y0vec);  free(z0vec);  free(x1vec);  free(y1vec);  free(z1vec);
    free(fx0vec); free(fy0vec); free(fz0vec); free(fx1vec); free(fy1vec); free(fz1vec);
    free(gammax); free(gammay); free(gammaz);
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
    if(__APPLE__) {
        // catlas_sset(int <numelem>, int <val>, float* <vec>, incx <stride>)
        catlas_sset(n, 0, vec, 1);
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

// EXTRA

float* setXLatticePosition(long n, float *x, float dx) {
    // fcc lattice construction
    long i, ai; long counter = 0; long unitcells = ceil(n/4);

    for(i = 0; i < unitcells; i++) {
        for(ai=0; ai<4; ai++) {
            // atom arrangement
            if(ai == 0) x[counter] = dx * i;
            if(ai == 1) x[counter] = dx * (i+0.5);
            if(ai == 2) x[counter] = dx * (i+0.5);
            if(ai == 3) x[counter] = dx * i;
            counter++;
            if(counter>n) break;
        }
        if(counter>n) break;
    }
    return x;
}

float* setYLatticePosition(long n, float *y, float dy) {
    // fcc lattice construction
    long j, aj; long counter = 0; long unitcells = ceil(n/4);

    for(j = 0; j < unitcells; j++) {
        for(aj=0; aj<4; aj++) {
            // atom arrangement
            if(aj == 0) y[counter] = dy * j;
            if(aj == 1) y[counter] = dy * (j+0.5);
            if(aj == 2) y[counter] = dy * j;
            if(aj == 3) y[counter] = dy * (j+0.5);
            counter++;
            if(counter>n) break;
        }
        if(counter>n) break;
    }
    return y;
}

float* setZLatticePosition(long n, float *z, float dz) {
    // fcc lattice construction
    long k, ak; long counter = 0; long unitcells = ceil(n/4);

    for(k = 0; k < unitcells; k++) {
        for(ak=0; ak<4; ak++) {
            // atom arrangement
            if(ak == 0) z[counter] = dz * k;
            if(ak == 1) z[counter] = dz * k;
            if(ak == 2) z[counter] = dz * (k+0.5);
            if(ak == 3) z[counter] = dz * (k+0.5);
            counter++;
            if(counter>n) break;
        }
        if(counter>n) break;
    }
    return z;
}

// IO
void PrintMatrix(const long n, const long m, float *vec) {
    // print matrix
    for(long i = 0; i < n*m; i++) {
        cout << std::scientific;
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

int Coords2XYZFile(long n, float *x, float *y, float *z, int index) {
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
        myfile << x[i] << "\t";
        myfile << y[i] << "\t";
        myfile << z[i] << "\n";
    }
    myfile.close();
    return 1;
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

float* ComputeGamma(long natoms, long ndim, float* xyzmat, float* xyzmat0, float* gammaxyz, float* fxyz, float* fxyz0) {
    // compute gamma
    float numerator_sum, denominator_sum, dxyz, df;
    long index; long i, j, k;
    for(i = 0; i < natoms; i++) {
        for(j = 0; j < ndim; j++) {
            numerator_sum = denominator_sum = 0;
            // inner product sum
            for(k = 0; k < natoms; k++) {
                dxyz = xyzmat[k*ndim + j] - xyzmat0[k*ndim + j];
                df = fxyz[k*ndim + j] - fxyz0[k*ndim + j];
                numerator_sum += dxyz * df; 
                denominator_sum += pow(df, 2);
            }
            index = i*ndim + j;
            gammaxyz[index] = numerator_sum / denominator_sum;
        }
    }
    return gammaxyz;
}

float* ComputeGammaVec(long natoms, float* x, float* x0, float* fx, float* fx0, float* gammax) {
    // compute gamma
    float numerator_sum, denominator_sum, dx, df;
    long index; long i, j, k;
    for(i = 0; i < natoms; i++) {
        numerator_sum = denominator_sum = 0;
        // inner product sum
        for(k = 0; k < natoms; k++) {
            dx = x[k] - x0[k];
            df = fx[k] - fx0[k];
            numerator_sum += dx * df;
            denominator_sum += pow(df, 2);
        }
        gammax[i] = numerator_sum / denominator_sum;
    }
    return gammax;
}

