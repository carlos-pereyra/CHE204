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
#include <algorithm>    // std::sort
#include <vector>       // std::vector

#if defined(__LP64__) /* In LP64 match sizes with the 32 bit ABI */
typedef int 		__CLPK_integer;
typedef int 		__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef int 		__CLPK_ftnlen;
#else
typedef long int 	__CLPK_integer;
typedef long int 	__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef long int 	__CLPK_ftnlen;
#endif

#define DBG 0
#define PNT 0

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
typedef __CLPK_doublereal doublereal;
typedef __CLPK_integer integer;
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

// INITIALIZATION
double* ReadInput(string filename, long nlines, double* mat);
double* ComputeKij(const long n, double* XYZi, double* Kij);

// VECTOR OPERATIONS
double* FormIdentityMatrix(long n, double* A);
double* VectorFill(const long n, double *v, double val);
float* VectorScalarProduct(long n, float a, float* Xi, float* Ri);
float* VectorAddition(long n, float* Xi, float* Yi, float* Ri);
float* VectorSubtraction(long n, float* Xi, float* Yi, float* Ri);
float  VectorDotProduct(long n, float* Xi, float* Yi);
float* VectorOuterProduct(long n, float* Xi, float* Yj, float* Rij);
int GetAbsMaxElementIndex(long n, long m, double* mat);

// MATRIX OPERATIONS
double* ClearMatrix(long n, long m, double* mat);
float* MatrixVectorProduct(long n, float* Aij, float* Xi, float* Ri);
float* MatrixMatrixProduct(long n, float* Aij, float* Bij, float* Cij);
float* MatrixCopy(long n, long m, float* x, float* y);
float* VectorExtraction(long n, long m, float*, long, float*);

// IO
void PrintMatrix(const long n, const long m, double *);
void PrintMatrixDiag(const long n, const long m, float *);
void Coords2XYZFile(long n, double* xyz, long index);
void Results2File(long iter, float e, float maxf);

// PHYSICS
double MorsePotential(double xi, double xj);
double ComputeTotalMorsePotential(long natoms, double *xyz);

// FORCES
double MorseForce(double xi, double xj, double rij);
double* ComputeMorseForces(long natoms, double *xyz, double *f);

// HESSIAN
double* ComputeHessian(long n, double* x, double* hij);

int main(int argc, char** argv) {
    // setup and initialization
    long natoms = 32; long ndim = 3; // make sure natoms matches the number of atoms in .xyz file
    long n = ndim * natoms;
    long numruns = 1100;
    double mem;

    /*
     * GEOMETRY OPTIMIZATION RELATED
     */
    //& coords
    double *XYZi = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;    // size (n by 1)
    double *SXYZi = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;   // size (n by 1)
    //& jacobian force matrix 
    double *Fj = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;      // size (n by 1)
    double *FjShift = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n; // size (n by 1)
    //& matrices
    double *Bjk = (double *) malloc(sizeof(double)*n*n); mem += sizeof(double)*n*n; // size (n by n)
    double *Cjk = (double *) malloc(sizeof(double)*n*n); mem += sizeof(double)*n*n; // size (n by n)
    double *Djk = (double *) malloc(sizeof(double)*n*n); mem += sizeof(double)*n*n; // size (n by n)
    double *Ejk = (double *) malloc(sizeof(double)*n*n); mem += sizeof(double)*n*n; // size (n by n)
    //& vectors
    double *Pk = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;      // size (n by 1)
    double *Sk = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;      // size (n by 1)
    double *Yk = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;      // size (n by 1)
    double *Dk = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;      // size (n by 1)
    //int *pivot = (int *) malloc(sizeof(int)*n);      // size (natom by ndim)
    
    /*
     * NORMAL MODES RELATED
     */
    double *Kij   = (double *) malloc(sizeof(double)*n*n); mem += sizeof(double)*n*n; // size (n by n)
    double *dXYZi = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;     // size (n by 1)
    double *Fi_ph = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;     // size (n by 1)
    double *Fi_mh = (double *) malloc(sizeof(double)*n); mem += sizeof(double)*n;     // size (n by 1)

    // read xyz positions and 
    ReadInput("coordinates/Configuration.xyz", natoms, XYZi);

    // initialize vectors
    FormIdentityMatrix(n, Bjk);
    VectorFill(n, Pk, 0);
    VectorFill(n, Yk, 0);
    VectorFill(n*n, Cjk, 0);
    VectorFill(n*n, Djk, 0);
    VectorFill(n*n, Kij, 0);

    // initialize scalars
    double ep, maxf;
    double a = 1; double b = 0; double alpha = 0.01;
    long stride = 1; long nshow = 5;

    /*
     * MINIMIZE POTENTIAL (GEOMETRY OPTIMIZATION)
     */
    for(long iter = 0; iter < numruns; iter++) {
        printf("\n=======================");
        printf("\n__ iteration = %ld __\n", iter);   
        printf("\n=======================\n");

        /*
         * step 0. clear all forces (n by 1)
         */
        ClearMatrix(n*n, 1, Cjk);
        ClearMatrix(n*n, 1, Djk);
        ClearMatrix(n*n, 1, Ejk);
        ClearMatrix(n, 1, Dk);
        ClearMatrix(n, 1, Fj);
        ClearMatrix(n, 1, FjShift);
        
        /*
         * step 1. compute fx, fy, fz vector (n by 1) 
         */
        ComputeMorseForces(natoms, XYZi, Fj);

        /* 
         * step 2. solve for Pk vector
         */
        int err; int nrhs = 1;
        cblas_dcopy(n, Fj, 1, Pk, 1);
        //dgesv_(&n, &nrhs, Bjk, &n, pivot, Pk, &n, &err); // find Pk in Bjk Pk = -\nabla f(Xi)
        if(err > 0) printf("factorization has been completed (U is exactly singular)\n");
        if(err < 0) printf("argument had an illegal value\n");
        if(err == 0) printf("successful"); 
        /*
         * step 3. solve for Sk = alpha Pk 
         */
        cblas_dcopy(n, Pk, 1, Sk, 1); // copy Pk to Sk
        cblas_dscal(n, alpha, Sk, 1); // Sk = a Sk

        /*
         * step 4. solve for Xk+1 = Xk + Sk evolved position
         */
        cblas_dcopy(n, Sk, 1, SXYZi, 1); // copy Sk to sXYZ
        cblas_daxpy(n, a, XYZi, stride, SXYZi, stride); // Xk+1 = Xk + Sk

        /*
         * step 5. solve for Yk = -FjShifted + Fj
         */
        ComputeMorseForces(natoms, SXYZi, FjShift);
        cblas_dcopy(n, Fj, 1, Yk, 1);
        cblas_daxpy(n, -a, FjShift, stride, Yk, stride); // Yk = -FjS + (Yk or Fj)

        /* step 6. 
            solve for Bjk' = Bjk + (Yk YkT) / (YkT Sk) - (Bjk Sk SkT BjkT) / (SjT Bjk Sj)
                      Bjk' = Bjk + Cjk / c             - (Bjk Sk Bjk Sk) / (SjT Dk)
                      Bjk' = Bjk + Cjk                 - (Bjk Sk Dk) / d
                      Bjk' = Bjk + Cjk                 - (Bjk Djk) / d
                      Bjk' = Bjk + Cjk                 - Ejk / d            
                      Bjk' = Bjk + Cjk                 - Ejk                */
        /*
         * step 6.1 solve for c (second term denominator) with dot product
         */
        double c = 0;
        c = cblas_ddot(n, Yk, stride, Sk, stride); // c = Yk * Sk
        
        /*
         * step 6.2 solve for Cjk (second term numerator)
         */
        ClearMatrix(n*n, 1, Cjk);
        cblas_dger(CblasRowMajor, n, n, a, Yk, stride, Yk, stride, Cjk, n); // Cjk = a Yk Yk^T + Cjk
        cblas_dscal(n, 1/c, Cjk, stride); // Cjk = Cjk / c

        /*
         * step 6.3 solve for Dk and then scalar d (third term denominator)
         */
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, a, Bjk, n, Sk, stride, b, Dk, stride); // Dk = a Bjk Sk + b Dk
        double d = 0;
        d = cblas_ddot(n, Sk, stride, Dk, stride); // d = Yk * Sk
        
        /*
         * step 6.4 solve for Ejk (third term numerator)
         */
        ClearMatrix(n, n, Djk);
        ClearMatrix(n, n, Ejk);
        cblas_dger(CblasRowMajor, n, n, a, Sk, stride, Dk, stride, Djk, n); // Djk = a Sk Dk + Djk
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, a, Bjk, n, Djk, n, b, Ejk, n); // Ejk = a Bjk Djk + b Ejk
        cblas_dscal(n*n, 1/d, Ejk, 1); // Ejk / d

        /*
         * step 7. Bj+1 = Bjk + a Cjk - a Djk 
         */
        cblas_daxpy(n*n, a, Cjk, stride, Bjk, stride);
        cblas_daxpy(n*n,-a, Ejk, stride, Bjk, stride);
        if(1) printf("\n__B_jk_FINAL__\n");
        if(1) PrintMatrix(nshow, 1, Bjk);

        /*
         * step 8. potential energy summary
         */
        ep = ComputeTotalMorsePotential(natoms, XYZi);
        printf("\n__EP__\n");
        printf("\nEP = %f\n", ep);

        /*
         * step 9. save shifted xyz to xyz
         */
        cblas_dcopy(n, SXYZi, 1, XYZi, 1);

        /* 
         * step 10. maximum force of all the components
         */
        int idmax = cblas_idamax(n, Fj, stride); //GetAbsMaxElementIndex(n, 1, Fj);
        maxf = Fj[idmax];
        printf("\n__MAX_F__");
        printf("\n F = %f\n", abs(maxf));

        Coords2XYZFile(natoms, XYZi, iter);
        if( !(n % 2) ) Results2File(iter, ep, abs(maxf));
        if(abs(maxf) < (1e-4) ) {
            break;
        }
    }

    /*
     * COMPUTE HESSIAN
     *          Kij = d''V / dxi dxj
     *
     * CENTRAL SECOND DERIVATIVE 
     *          Kij = (-Fi(xi+h) - -Fi(xi-h)) / dx
     * 
     */
    ComputeHessian(n, XYZi, Kij);
    if(1) printf("\n__H_ij__\n");
    if(1) PrintMatrix(n*n, 1, Kij);
    
    /*
     * FREE MEMORY
     */
    printf("\nfree %lf (kB)\n", mem / 1000.);

    free(XYZi); free(SXYZi);    // GEOMETRY OPTIMIZATION
    free(Fj); free(FjShift);
    free(Bjk); free(Cjk); free(Djk); free(Ejk);
    free(Pk);  free(Sk); free(Yk); free(Dk);
    
    free(Kij);                  // NORMAL MODES 
    free(dXYZi); 
    free(Fi_ph); free(Fi_mh);
}

// INITIALIZATION

double* ReadInput(string filename, long nlines, double *mat) {
    // read in positions from .xyz file
    ifstream infile(filename);
    string elem, line;
    double x, y, z;
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

        mat[3*counter + 0] = x;
        mat[3*counter + 1] = y;
        mat[3*counter + 2] = z;

        counter++;
    }
    infile.close();
    return mat;
}

double* ComputeKij(const long n, double* XYZi, double* Kij) {
    /*
     * Compute the constant force matrix Kij
     */
    long i, j, k; double xi, xj; double h = 0.01; double Vxij, Vxij_ph, Vxij_mh;
    //std::vector<long> myvector;
    long ndim = 3;
    long natom = n / (double) ndim;

    // particle i
    
    for(i = 0; i < natom; i++) {
        
        // particle j
        
        for(j = 0; j < natom; j++) {
            
            // 3 dimensions
        
            for(k = 0; k < ndim; k++) {
                xi = XYZi[i*3 + k];
                xj = XYZi[j*3 + k];
                printf("ComputeKij = %ld\n", i*n + j*3 + k );
                Vxij    = MorsePotential(xi, xj);
                Vxij_ph = MorsePotential(xi + h, xj);
                Vxij_mh = MorsePotential(xi - h, xj);
                //vector.push_back(i*(j + k) );
                //Kij[i*n + j]
            }
        }
    }
    return Kij;
}

// MATRIX OPERATIONS

double* ClearMatrix(long n, long m, double *vec) {
    // clear vector that is a real n by m matrix
    if(DBG) printf("\nClearMatrix()\n");
    if(__APPLE__) {
        long stride = 1;
        catlas_dset(n*m, 0, vec, stride);
    } else {
        for(long i = 0; i < n; i++) {
            for(long j = 0; j < m; j++) {
                vec[i*n + j] = 0;
            }
        }
    }
    return vec;
}

int GetAbsMaxElementIndex(long n, long m, double *vec) {
    // get abs(max) element in vector that is a real n by m matrix
    // skip is the stride constant
    int index; float max;
    if(__APPLE__) {
        long stride = 1;
        index = cblas_idamax(n*m, vec, stride);
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
double* VectorFill(const long n, double *v, double val) {
    // fill vector with single value
    // return filled vector
    //if(__APPLE__) {
    //    catlas_sset(n, val, v, 1);
    //} else {
        for(long i=0; i<n; i++) {
            v[i] = val;
        }
    //}
    return v;
}

double* FormIdentityMatrix(long n, double* A) {
    // form an identity matrix
    long diag, index;
    for(long i = 0; i < n; i++) {
        diag = i * (n + 1);
        for(long j = 0; j < n; j++) {
            index = i*n + j;
            if(index == diag) A[index] = 1;
            else A[index] = 0;
        }
    }
    return A;
}

float* VectorScalarProduct(long n, float a, float* Xi, float* Ri) {
    /* Scale a vector Xi by scalar a and return result vector Ri */
    if(__APPLE__) {
        cblas_scopy(n, Xi, 1, Ri, 1);
        cblas_sscal(n, a, Ri, 1);
    } else {
        for(long i = 0; i < n; i++) {
            Ri[i] = a * Xi[i];
        }
    }
    return Ri;
}

float* VectorAddition(long n, float* Xi, float* Yi, float* Ri) {
    /* sum vectors Xi and Yi into Ri vector. be aware Yi copies into Ri. */
    float scalar = 1; float inc = 1;
    if(__APPLE__) {
        cblas_scopy(n, Yi, 1, Ri, 1);
        cblas_saxpy(n, scalar, Xi, inc, Ri, inc);
    } else {
        for(long i = 0; i < n; i++) {
            Ri[i] = Yi[i] + Xi[i];
        }
    }
    return Ri;
}

float* VectorSubtraction(long n, float* Xi, float* Yi, float* Ri) {
    /* subtract Yi elements by Xi into Ri vector. be aware Yi copies into Ri.
                Y[i] = alpha * X[i] + Y[i]

       saxpy	(	integer 	N,
                    real 	SA,
                    real, dimension(*) 	SX,
                    integer 	INCX,
                    real, dimension(*) 	SY,
                    integer 	INCY)	*/
    float scalar = -1; float stride = 1;
    if(__APPLE__) {
        cblas_scopy(n, Yi, 1, Ri, 1);
        cblas_saxpy(n, scalar, Xi, stride, Ri, stride);
    } else {
        for(long i = 0; i < n; i++) {
            Ri[i] = Yi[i] - Xi[i];
        }
    }
    return Ri;
}   

float VectorDotProduct(long n, float* Xi, float* Yi) {
    /* vector Xi and Yi dot product returns scalar d */
    float sum; float result;
    if(__APPLE__) {
        long stride = 1;
        result = cblas_sdot(n, Xi, stride, Yi, stride);
    } else {
        sum = 0;
        for(long i = 0; i < n; i++) {
            sum += Xi[i] * Yi[i];
        }
        result = sum;
    }
    return result;
}

float* VectorOuterProduct(long n, float* Xi, float* Yj, float* Rij) {
    /* outer product between Xi and Yj return Rij matrix
        cblas_sger solves,
                            Rij = alpha Xi Yj + Rij 
        sger(integer M,
             integer N,
             real 	ALPHA,
             real, dimension(*) 	X,
             integer 	INCX,
             real, dimension(*) 	Y,
             integer 	INCY,
             real, dimension(lda,*) 	A,
             integer 	LDA)	                        */
    long i, j;
    if(__APPLE__) {
        float alpha = 1; long stride = 1; 
        cblas_sger(CblasRowMajor, n, n, alpha, Xi, stride, Yj, stride, Rij, n);
    } else {
        for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                Rij[i*n + j] = Xi[i] * Yj[j];
            }
        }
    }
    return Rij;
}

float* MatrixVectorProduct(long n, float* Aij, float* Xi, float* Ri) {
    /* matrix Aij and vector Xi product result is Ri 
       cblas_sgemv solves,
                            Yi = a Aij Xi + b Yi 

       sgemv	(	character 	TRANS,
                    integer 	M,
                    integer 	N,
                    real 	ALPHA,
                    real, dimension(lda,*) 	A,
                    integer 	LDA,
                    real, dimension(*) 	X,
                    integer 	INCX,
                    real 	BETA,
                    real, dimension(*) 	Y,
                    integer 	INCY)	*/
    float sum; long i, j;
    if(__APPLE__) {
        float alpha = 1; long stride = 1; float beta = 0;
        cblas_scopy(n, Xi, 1, Ri, 1);
        cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, alpha, Aij, n, Xi, stride, beta, Ri, stride);
    } else {
        for(i = 0; i < n; i++) {
            sum = 0;
            for(j = 0; j < n; j++) {
                sum += Aij[i*n + j] * Xi[j];
            }
            Ri[i] = sum;
        }
    }
    return Ri;
}

float* MatrixMatrixProduct(long n, float* Aij, float* Bij, float* Cij) {
    /* matrix matrix product 
       Cij = a Aij Bij + b Cij

       sgemm	(	character 	TRANSA,
                    character 	TRANSB,
                    integer 	M,
                    integer 	N,
                    integer 	K,
                    real 	ALPHA,
                    real, dimension(lda,*) 	A,
                    integer 	LDA,
                    real, dimension(ldb,*) 	B,
                    integer 	LDB,
                    real 	BETA,
                    real, dimension(ldc,*) 	C,
                    integer 	LDC )	 */
    long i, j; float sum;
    if(__APPLE__) {
        float alpha = 1; long stride = 1; float beta = 0;
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, Aij, n, Bij, n, beta, Cij, n);
    } else {
        /*for(i = 0; i < n; i++) {
            sum = 0;
            for(j = 0; j < n; j++) {
                sum += Aij[i*n + j] * Bij[j*n + i];
                Cij[i*n + j] = sum;
            }
        }*/
    }
    return Cij;
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
void PrintMatrix(const long n, const long m, double *vec) {
    // print matrix
    for(long i = 0; i < n*m; i++) {
        //cout << std::scientific;
        cout.precision(4);
        cout << setw(8) << i << "\t "<< vec[i] << "\t\n";
        //if( !((i+1)%m) ) { cout << "\n"; }
    }
}

void Coords2XYZFile(long n, double *xyz, long index) {
    // write coordinates to xyz file
    if(DBG) printf("\nCoords2XYZFile()\n");
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

/*
 * POTENTIAL CALCULATION
 */

double MorsePotential(double xi, double xj) {
    // separation distance between i'th and j'th particle = 'float rij'
    // magnitude of separation between i'th and j'th particle = 'float drij'
    double d = 1; double a = 1; double xo = 3.5;
    double xij = xi - xj;
    double ep = d * pow((1 - exp(-a*(xij-xo))), 2);
    return ep;
}

double ComputeTotalMorsePotential(long natoms, double *xyz) { 
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
            ep_sum += ep_value;
        }
        global_ep_sum += ep_sum;
    }
    return global_ep_sum;
}

/*
 * FORCE CALCULATIONS
 */

double MorseForce(double xi, double xj, double rij) {
    //
    if(DBG) printf("\nMorseForce()\n");
    double d = 1; double a = 1; double ro = 3.5;
    double xij = xi - xj;
    double f = 2*d*a * ( exp(-2*a*(rij-ro)) - exp(-a*(rij-ro)) ) * xij / rij;
    return f;
}

double* ComputeMorseForces(long natoms, double *xyz, double *f) {
    if(DBG) printf("\nComputeMorseForces()\n");
    float d = 1; float a = 1; float ro = 3.5;
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


/*
 * HESSIAN
 */

double* ComputeHessian(long n, double* x, double* hij) {
    // compute 3N x 3N real hessian matrix
    // finite central difference method of second order derivative requires differential distance = 'float d'
    // xyz positions = 'float* xyz'
    // 3n x 3n hessian matrix = 'float* hij'
    long i, j, k; long count = 0;
    double xi, xj;
    double dx = 0.01;
    double xij, yij, zij;
    double rij_plus;
    double rij_minus;
    double f1, f2;
    /*
     * F(x + dx) = f1 (force)
     * F(x - dx) = f2 (force)
     */
    long natom = n / 3;
    for(i = 0; i < n; i++) {
        for(j = 0; j < natom; j++) {
            double xij       = x[3*i + 0] - x[3*j + 0];
            double yij       = x[3*i + 1] - x[3*j + 1];
            double zij       = x[3*i + 2] - x[3*j + 2];
            double rij_plus  = sqrt(pow(xij+dx, 2) + pow(yij+dx, 2) + pow(zij+dx, 2));
            double rij_minus = sqrt(pow(xij-dx, 2) + pow(yij-dx, 2) + pow(zij-dx, 2));
            for(k = 0; k < 3; k++) {
                xi = x[3*i + k];
                xj = x[3*j + k];
                f1 = MorseForce(xi+dx, xj, rij_plus);
                f2 = MorseForce(xi-dx, xj, rij_minus);
                hij[i*n + j*3 + k] = ( -1*f1 - -1*f2 ) / dx;
                //printf("i = %ld j = %ld k = %ld xi = %lf xj = %lf f =%lf\n", i, j, k, xi, xj, hij[i*n + j*3 + k]);
                //count++;
            }
        }
    }
    return hij;
}
