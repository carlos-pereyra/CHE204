#include <cstdlib>    // c++
#include <cstdio>     // c++
#include <omp.h>
#include <vector>     // c++
#include <algorithm>  // c++
#include <iostream>   // c++
#include <iomanip>    // c++
#include <fstream>    // c++
#include <math.h>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif

// (on OS X) compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -framework Accelerate -o test main.c
//
//      ./test
//

#define ranfloat(w) ( w * ((rand() / (double) RAND_MAX)  - 0.5) )
#define DBG 0
using namespace std;

float* readInput(string filename, long nlines, float *mat);
float* ClearMatrix(long n, long m, float* mat);
int GetAbsMaxElementIndex(long n, long m, float* mat, long skip);
// EXTRA
float* setSMatrix(long n, long m, float sxy, float syz, float szy); // ignore
float* setXLatticePosition(long n, float *x, float dx);             // ignore
float* setYLatticePosition(long n, float *y, float dy);             // ignore
float* setZLatticePosition(long n, float *z, float dz);             // ignore

// IO
void PrintMatrix(const long n, const long m, float *, long skip);
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

int main(int argc, char** argv) {
    long natoms = 3; long ndim = 3; // make sure natoms matches the number of atoms in .xyz file
    // initialize
    float *xyzmat = (float *) malloc(sizeof(float)*natoms*ndim);    // size (natom by 3)
    
    float *dxmat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *dymat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *dzmat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *drmat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
   
    float *fxmat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *fymat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    float *fzmat = (float *) malloc(sizeof(float)*natoms*natoms);   // size (natom by natom)
    // read xyz positions
    readInput("coordinates/Test.xyz", natoms, xyzmat);
    printf("XYZ__\n"); PrintMatrix(natoms, 3, xyzmat, 1);
    //readInput("coordinates/Configuration.xyz", natoms, xyzmat);
    // clear all forces
    ClearMatrix(natoms, natoms, fxmat);
    ClearMatrix(natoms, natoms, fymat);
    ClearMatrix(natoms, natoms, fzmat);
    // compute dx, dy, dz, dr
    ComputeDX(natoms, xyzmat, dxmat);
    ComputeDY(natoms, xyzmat, dymat);
    ComputeDZ(natoms, xyzmat, dzmat);
    ComputeDR(natoms, xyzmat, drmat);
    // compute fx, fy, fz jacobian matrices (diagonal elements are atom forces)
    ComputeMorseForceX(natoms, dxmat, drmat, fxmat);
    ComputeMorseForceY(natoms, dymat, drmat, fymat);
    ComputeMorseForceZ(natoms, dzmat, drmat, fzmat);
    printf("\nFX__\n"); PrintMatrix(natoms, natoms, fxmat, natoms+1);
    printf("\nFY__\n"); PrintMatrix(natoms, natoms, fymat, natoms+1);
    printf("\nFZ__\n"); PrintMatrix(natoms, natoms, fzmat, natoms+1);
    int idx1 = GetAbsMaxElementIndex(natoms, natoms, fxmat, natoms+1);
    printf("idx1 = %d\n", idx1);

    // clean
    free(xyzmat); free(dxmat); free(dymat); free(dzmat); free(drmat); 
    free(fxmat); free(fymat); free(fzmat);
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

int GetAbsMaxElementIndex(long n, long m, float *vec, long skip) {
    // get abs(max) element in vector that is a real n by m matrix
    // skip is the stride constant
    int index; float max;
    if(__APPLE__) {
        printf("cblas_isamax works!\n");
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
void PrintMatrix(const long n, const long m, float *vec, const long skip) {
    // print matrix
    for(long i = 0, j = 0; j < n*m; i++, j=i*skip) {
        //if( i*skip > n*m ) break;
        //printf(" i = %ld ", i);
        cout << std::scientific;
        cout.precision(2);
        cout << setw(8) << vec[i*skip] << "\t";
        if( !((i+1)%n) ) { cout << "\n"; }
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
