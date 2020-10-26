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

// VECTOR OPERATIONS
float* VectorFill(const long n, float *v, float val);

// GDA RELATED
float* SetInitialGuessPositions(long n, long m, MTRand *, float* xyz, float* newxyz);
float* ComputeGammaX(long n, long m, float* xyzmat, float* xyzmat0, float* gamma);
float* MatrixCopy(long n, long m, float* x, float* y);
float* VectorExtraction(long n, long m, float*, long, float*);

// IO
void PrintMatrix(const long n, const long m, float *);
void PrintMatrixDiag(const long n, const long m, float *);
void Coords2XYZFile(long n, float* xyz, long index);
void Results2File(long iter, float e, float maxf);

// PHYSICS
float* ComputeHessian(long natom, long ndim, float d, float* xyz, float* hij);

float ComputeAMorsePotential(float rij);
float ComputeTotalMorsePotential(long natoms, float *xyz);

//      a. FORCES
float* ComputeMorseForces(long natoms, float *xyz, float *f);
float* ComputeMorseForceX(long n, float *xyz, float *fx);
float* ComputeMorseForceY(long n, float *xyz, float *fy);
float* ComputeMorseForceZ(long n, float *xyz, float *fz);

//      b. POSITION
float* ComputeNewX(long n, float* xyzmat, float* g, float* fxmat);
float* ComputeNewY(long n, float* xyzmat, float* g, float* fymat);
float* ComputeNewZ(long n, float* xyzmat, float* g, float* fzmat);

float* ComputeNewGamma(long n, float* xyz, float* xyzold, float* g, float* fx, float* fxo, float* fy, float* fyo, float* fz, float* fzo);

int main(int argc, char** argv) {
    // setup and initialization
    long natoms = 32; long ndim = 3;    // make sure natoms matches the 
                                       // number of atoms in .xyz file
    long numruns = 1;

    //& coords
    float *xyzmat_old = (float *) malloc(sizeof(float)*natoms*ndim);    // size (natom by ndim)
    float *xyzmat = (float *) malloc(sizeof(float)*natoms*ndim);        // size (natom by ndim)
    //& pair distance
    float *dxvec = (float *) malloc(sizeof(float)*natoms);      // size (natom by natom)
    float *dyvec = (float *) malloc(sizeof(float)*natoms);      // size (natom by natom)
    float *dzvec = (float *) malloc(sizeof(float)*natoms);      // size (natom by natom)
    float *drvec = (float *) malloc(sizeof(float)*natoms);      // size (natom by natom)
    //& jacobian force matrix 
    float *Fvec_old = (float *) malloc(sizeof(float)*ndim*natoms);      // size (natom by 1)
    float *fxvec_old = (float *) malloc(sizeof(float)*natoms);  // size (natom by 1)
    float *fyvec_old = (float *) malloc(sizeof(float)*natoms);  // size (natom by 1)
    float *fzvec_old = (float *) malloc(sizeof(float)*natoms);  // size (natom by 1)
    float *Fvec = (float *) malloc(sizeof(float)*ndim*natoms);  // size (natom by 1)
    float *fxvec = (float *) malloc(sizeof(float)*natoms);      // size (natom by 1)
    float *fyvec = (float *) malloc(sizeof(float)*natoms);      // size (natom by 1)
    float *fzvec = (float *) malloc(sizeof(float)*natoms);      // size (natom by 1)
    float *maxfvec = (float *) malloc(sizeof(float)*ndim);      // size (ndim by 1)
    //& hessian matrix
    float *hij = (float *) malloc(sizeof(float)*ndim*natoms*ndim*natoms);// size (natom by 1)
    //& gamma vector
    float *gamma = (float *) malloc(sizeof(float)*natoms);      // size (natom by 1)
    float *g = (float *) malloc(sizeof(float)*natoms);          // size (natom by 1)
    // random number generator
    unsigned long int now = static_cast<unsigned long int>(time(NULL));
    MTRand *mtrand = new MTRand(now);

    // read xyz positions and 

    //readInput("coordinates/Test.xyz", natoms, xyzmat);
    readInput("coordinates/Configuration.xyz", natoms, xyzmat);

    // initial guess gamma values 
    VectorFill(natoms, gamma, 0.01);
    VectorFill(natoms, g, 0.01);
    float ep, maxf;
    
    for(long n=0; n<numruns; n++) {
        printf("\n=======================");
        printf("\n__ iteration = %ld __\n", n);   
        printf("\n=======================\n");
        MatrixCopy(natoms, ndim,   xyzmat, xyzmat_old);
        MatrixCopy(natoms, ndim, Fvec,  Fvec_old);

        // clear all forces (n x 1)
        ClearMatrix(natoms, 1, Fvec);
        
        // compute total system energy
        ep = ComputeTotalMorsePotential(natoms, xyzmat);
        
        // compute fx, fy, fz vector (n x 1)
        ComputeMorseForces(natoms, xyzmat, Fvec);

        printf("\n__EP__\n");
        printf("\n EP = %f\n", ep);
        printf("\n__F_VECTOR__\n");
        PrintMatrix(natoms, ndim, Fvec);

        // compute new positions (returns new xyzmat)
        /*ComputeNewX(natoms, xyzmat, gamma, fxvec); 
        ComputeNewY(natoms, xyzmat, gamma, fyvec);
        ComputeNewZ(natoms, xyzmat, gamma, fzvec);
        */
        
        // maximum force of all the components
        int idMax = GetAbsMaxElementIndex(ndim, 1, Fvec, 1);
        maxf = Fvec[idMax];

        // show absolute maximum element
        //cout << std::scientific;
        cout.precision(2);
        cout << "================================";
        cout << "\nAbs. Max force component is = ";
        cout << Fvec[idMax] << "\n";

        Coords2XYZFile(natoms, xyzmat, n);
        if( !(n % 2) ) Results2File(n, ep, maxf);

        if(abs(maxf) < (1e-4) ) {
            break;
        }
    }

    // free memory 
    free(xyzmat); free(xyzmat_old);
    free(dxvec);    free(dyvec);     free(dzvec);     free(drvec); 
    free(Fvec);     free(fxvec);     free(fyvec);     free(fzvec); 
    free(Fvec_old); free(fxvec_old); free(fyvec_old); free(fzvec_old);
    free(hij);
    free(gamma);
    free(g);
    free(maxfvec);
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

// VECTOR OPERATIONS
float* VectorFill(const long n, float *v, float val) {
    // fill vector with single value
    // return filled vector
    if(__APPLE__) {
        catlas_sset(n, val, v, 1);
    } else {
        for(int i=0; i<=n; i++) {
            v[i] = val;
        }
    }
    return v;
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
        cout << setw(8) << vec[i] << "\t\n";
        //if( !((i+1)%m) ) { cout << "\n"; }
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

void Results2File(long iter, float e, float maxf) {
    // write energy and maximum force to file
    string filename = "data/output.dat";
    ofstream myfile;
    long i;
    // open new file
    if (iter == 0) {
        myfile.open(filename, std::fstream::out);
        myfile << left << setw(20) << "# iter";
        myfile << left << setw(20) << "e";
        myfile << "maxf" << "\n";
    }
    // append to file
    else myfile.open(filename, std::fstream::app);

    // write results to file
    myfile << std::scientific;
    myfile << left << setw(20) << iter;
    myfile << left << setw(20) << e;
    myfile << maxf << "\n";
    myfile.close();
}

// SEPARATION DISTANCES

float* ComputeHessian(long natom, long ndim, float d, float* xyz, float* hij) {
    // compute 3N x 3N real hessian matrix
    // finite central difference method of second order derivative requires differential distance = 'float d'
    // xyz positions = 'float* xyz'
    // 3n x 3n hessian matrix = 'float* hij'
    long i, j, k;
    for(i = 0; i < natom; i++) {
        for(j = 0; j < natom * ndim; j++) {
            //float dXij   = xyz[3*i + 0] - xyz[3*j + 0];
            //float dYij   = xyz[3*i + 1] - xyz[3*j + 1];
            //float dZij   = xyz[3*i + 2] - xyz[3*j + 2];
            //float dRij   = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
            //float dRijpd = sqrt(pow(dx+d, 2) + pow(dy+d, 2) + pow(dz+d, 2));
            //float dRijmd = sqrt(pow(dx-d, 2) + pow(dy-d, 2) + pow(dz-d, 2));
            for(k = 0; k < 3; k++) {
                float Rij = xyz[3*i + k] - xyz[3*j + k];
                float Eppd = ComputeAMorsePotential(Rij + d);
                float Epmd = ComputeAMorsePotential(Rij - d);
                float Ep = ComputeAMorsePotential(Rij);
                hij[3*natom*i + 3*j + k] = ( Eppd + Epmd - 2*Ep ) / pow(d, 2);
            }
        }
    }
    return hij;
}

// POTENTIAL CALCULATION

float ComputeAMorsePotential(float rij) {
    // separation distance between i'th and j'th particle = 'float rij'
    // magnitude of separation between i'th and j'th particle = 'float drij'
    float d = 1; float a = 1; float ro = 3.5;
    float ep = d * pow((1 - exp(-a*(rij-ro))), 2);
    return ep;
}

float ComputeTotalMorsePotential(long natoms, float *xyz) { //, float *ep) {
    // this function computes the morse potential. 
    // potential energy vector (spans natoms) = 'float* ep' 
    //
    // return total potential energy.
    float d = 1; float a = 1; float ro = 3.5;
    float ep_value; float ep_sum; float global_ep_sum = 0;
    float dx, dy, dz, dr;
    long i, j;
    for(i = 0; i < natoms; i++) {
        ep_sum = 0;
        for(j = i + 1; j < natoms; j++) {
            dx = xyz[i*3 + 0] - xyz[j*3 + 0];
            dy = xyz[i*3 + 1] - xyz[j*3 + 1];
            dz = xyz[i*3 + 2] - xyz[j*3 + 2];
            dr = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
            ep_value = d * pow((1 - exp(-a*(dr-ro))), 2);
            //ep[j] += ep_value;
            ep_sum += ep_value;
        }
        global_ep_sum += ep_sum;
        //ep[i] = ep_sum;
    }
    return global_ep_sum;
}

// FORCE CALCULATIONS

float* ComputeMorseForces(long natoms, float *xyz, float *f) {
    if(DBG) printf("\nComputeMorseForces()\n");
    float d = 1; float a = 1; float ro = 3.5;
    //float d = 2; float a = 1; float ro = 1.2;
    float force; float force_sum;
    float dXij, dYij, dZij, dRij, Rij;
    long i, j, k;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        force_sum = 0;
        for(j = i + 1; j < natoms; j++) {
            dXij = xyz[i*3 + 0] - xyz[j*3 + 0];
            dYij = xyz[i*3 + 1] - xyz[j*3 + 1];
            dZij = xyz[i*3 + 2] - xyz[j*3 + 2];
            dRij = sqrt(pow(dXij, 2) + pow(dYij, 2) + pow(dZij, 2));
            for(k = 0; k < 3; k++) {
                Rij = xyz[i*3 + k] - xyz[j*3 + k];
                force = 2*d*a * ( exp(-2*a*(dRij-ro)) - exp(-a*(dRij-ro)) ) * Rij / dRij;
                // 
                //printf("particle %ld interacts w %ld in %ld dim, force = %f\n", i, j, k, force);
                f[3*i + k] += force;
                f[3*j + k] -= force;
            }
        }
    }
    return f;
}


float* ComputeMorseForceX(long natoms, float *xyz, float *fx) {
    // provide natoms for the shape of the following matrices:
    //      - dx matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fx matrix size (natoms x natoms)
    // return fx
    if(DBG) printf("\nComputeMorseForceX()\n");
    float d = 1; float a = 1; float ro = 3.5;
    float val = 0; float force_sum;
    float dx, dy, dz, dr;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        force_sum = 0;
        for(j = i + 1; j < natoms; j++) {
            dx = xyz[i*3 + 0] - xyz[j*3 + 0];
            dy = xyz[i*3 + 1] - xyz[j*3 + 1];
            dz = xyz[i*3 + 2] - xyz[j*3 + 2];
            dr = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
            val = 2*d*a * ( exp(-2*a*(dr-ro)) - exp(-a*(dr-ro)) ) * dx / dr;
            fx[j] -= val;
            force_sum += val;
        }
        // particle force elements
        fx[i] += force_sum;
    }
    return fx;
}

float* ComputeMorseForceY(long natoms, float *xyz, float *fy) {
    // provide natoms for the shape of the following matrices:
    //      - dy matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fx matrix size (natoms x natoms)
    // return fy matrix
    if(DBG) printf("\nComputeMorseForceY()\n");
    float d = 1; float a = 1; float ro = 3.5; 
    float val = 0; float force_sum;
    float dx, dy, dz, dr;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        force_sum = 0;
        for(j = i + 1; j < natoms; j++) {
            dx = xyz[i*3 + 0] - xyz[j*3 + 0];
            dy = xyz[i*3 + 1] - xyz[j*3 + 1];
            dz = xyz[i*3 + 2] - xyz[j*3 + 2];
            dr = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
            val = 2*d*a * ( exp(-2*a*(dr-ro)) - exp(-a*(dr-ro)) ) * dy / dr;
            fy[j] -= val;
            force_sum += val;
        }
        fy[i] += force_sum;
    }
    return fy;
}

float* ComputeMorseForceZ(long natoms, float *xyz, float *fz) {
    // provide natoms for the shape of the following matrices:
    //      - dz matrix size (natoms x natoms)
    //      - dr matrix size (natoms x natoms)
    //      - fz matrix size (natoms x natoms)
    // return fz matrix
    if(DBG) printf("\nComputeMorseForceZ()\n");
    float d = 1; float a = 1; float ro = 3.5;
    float val = 0; float force_sum;
    float dx, dy, dz, dr;
    long i, j;
    // upper diagonal forces in matrix
    for(i = 0; i < natoms; i++) {
        force_sum = 0;
        for(j = i + 1; j < natoms; j++) {
            dx = xyz[i*3 + 0] - xyz[j*3 + 0];
            dy = xyz[i*3 + 1] - xyz[j*3 + 1];
            dz = xyz[i*3 + 2] - xyz[j*3 + 2];
            dr = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
            val = 2*d*a * ( exp(-2*a*(dr-ro)) - exp(-a*(dr-ro)) ) * dz / dr;
            fz[j] -= val;
            force_sum += val;
        }
        fz[i] += force_sum;
    }
    return fz;
}

// NEW POSITIONS

float* ComputeNewX(long n, float* xyzmat, float* g, float* fxvec) {
    // compute new x coorindate with new gamma and force computations
    //      eqn: compute x^(k+1)[i] = x[i] + gammax[i]*fx[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewX()\n");
    
    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 0] = xyzmat[i*ndim + 0] + g[i] * fxvec[i];
    }
    return xyzmat;
}

float* ComputeNewY(long n, float* xyzmat, float* g, float* fyvec) {
    // compute new y coorindate with new gamma and force computations
    //      eqn: compute y^(k+1)[i] = y[i] + gy * fy[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewY()\n");
    
    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 1] = xyzmat[i*ndim + 1] + g[i] * fyvec[i];
    }
    return xyzmat;
}

float* ComputeNewZ(long n, float* xyzmat, float* g, float* fzvec) {
    // compute new z coorindate with new gamma and force computations
    //      eqn: compute z^(k+1)[i] = z[i] + gz * fz[i]
    // return new xyz (n x 3) coordinates
    if(DBG) printf("\nComputeNewZ()\n");

    long ndim = 3;
    for(long i = 0; i<n; i++) {
        xyzmat[i*ndim + 2] = xyzmat[i*ndim + 2] + g[i] * fzvec[i];
    }
    return xyzmat;
}

// NEW GAMMA

float* ComputeNewGamma(long n, float* xyz, float* xyzold, float* g, float* fx, float* fxo, float* fy, float* fyo, float* fz, float* fzo) {
    // compute new gamma for x direction
    //  expects xyz (natoms x 3) matrix
    //  expects xyzold (natoms x 3) matrix
    //  expects g (natoms x 1) vec
    //  expects fx (natoms x 1) vector
    //  expects fxold (natoms x 1) vector
    //  expects fy (natoms x 1) vector
    //  expects fyold (natoms x 1) vector
    //  expects fz (natoms x 1) vector
    //  expects fzold (natoms x 1) vector
    if(DBG) printf("\nComputeNewGamma()\n");
    long i;
    float dx, dy, dz, dfx, dfy, dfz;
    float numerator = 0;
    float denominator = 0;
    for(i = 0; i < n; i++) {
        dx = xyz[ i*3 + 0 ] - xyzold[ i*3 + 0 ];
        dy = xyz[ i*3 + 1 ] - xyzold[ i*3 + 1 ];
        dz = xyz[ i*3 + 2 ] - xyzold[ i*3 + 2 ];
        dfx = fx[i] - fxo[i];
        dfy = fy[i] - fyo[i];
        dfz = fz[i] - fzo[i];
        numerator = abs(dx * dfx + dy * dfy + dz * dfz);
        denominator = pow(dfx, 2) + pow(dfy, 2) + pow(dfz, 2);
        if((numerator == 0) || (denominator == 0)) {
             g[i] = 0;
             continue;
        }
        g[i] = numerator / denominator;
        //printf("g[%ld] = %f\n", i, g[i]);
    }
    // return vector
    return g;
}

