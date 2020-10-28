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

float* readInput(string filename, long nlines, float *mat);
float* ClearMatrix(long n, long m, float* mat);
float* SetMatrixElement(long n, long m, float *mat, long index, long value);
int GetAbsMaxElementIndex(long n, long m, float* mat, long skip);

// VECTOR OPERATIONS
float* VectorFill(const long n, float *v, float val);

// LU Decomposition
float* FormIdentityMatrix(long n, float* A);
float* LUDecomposition(long n, float* bjk);
float* SolvePk(long n, float* Bjk, float* Pk, float* Fj);

float* VectorScalarProduct(long n, float a, float* Xi, float* Ri);
float* VectorAddition(long n, float* Xi, float* Yi, float* Ri);
float* VectorSubtraction(long n, float* Xi, float* Yi, float* Ri);
float  VectorDotProduct(long n, float* Xi, float* Yi);
float* VectorOuterProduct(long n, float* Xi, float* Yj, float* Rij);

float* MatrixVectorProduct(long n, float* Aij, float* Xi, float* Ri);
float* MatrixMatrixProduct(long n, float* Aij, float* Bij, float* Cij);

// GDA RELATED
float* SetInitialGuessPositions(long n, long m, MTRand *, float* xyz, float* newxyz);
float* MatrixCopy(long n, long m, float* x, float* y);
float* VectorExtraction(long n, long m, float*, long, float*);

// IO
void PrintMatrix(const long n, const long m, float *);
void PrintMatrixDiag(const long n, const long m, float *);
void Coords2XYZFile(long n, float* xyz, long index);
void Results2File(long iter, float e, float maxf);

// PHYSICS
float* MorsePotential(long s, long e, float *r, float *ep);

float* ComputeHessian(long natom, long ndim, float d, float* xyz, float* hij);

float ComputeAMorsePotential(float rij);
float ComputeTotalMorsePotential(long natoms, float *xyz);

//      a. FORCES
float* ComputeMorseForces(long natoms, float *xyz, float *f);

//      b. POSITION
float* ComputeNewPositions(long natom_x_ndim, float* xyzvec, float* g, float* fvec);


int main(int argc, char** argv) {
    // setup and initialization
    long natoms = 32; long ndim = 3; // make sure natoms matches the number of atoms in .xyz file
    int n = ndim * natoms;
    long numruns = 9;

    double *vec = (double*) malloc(n*sizeof(double));
    double *mat = (double*) malloc(n*n*sizeof(double));
    //& coords
    float *sXYZ = (float *) malloc(sizeof(float)*n);    // size (natom by ndim)
    float *XYZ = (float *) malloc(sizeof(float)*n);     // size (natom by ndim)
    //& jacobian force matrix 
    float *Fj = (float *) malloc(sizeof(float)*n);      // size (natom by ndim)
    float *FjShift = (float *) malloc(sizeof(float)*n); // size (natom by ndim)
    //& matrices
    float *Bjk = (float *) malloc(sizeof(float)*n*n);   // size (natom*ndim by natom*ndim)
    float *Cjk = (float *) malloc(sizeof(float)*n*n);   // size (natom*ndim by natom*ndim)
    float *Djk = (float *) malloc(sizeof(float)*n*n);   // size (natom*ndim by natom*ndim)
    float *Ejk = (float *) malloc(sizeof(float)*n*n);   // size (natom*ndim by natom*ndim)
    //& vectors
    float *Pk = (float *) malloc(sizeof(float)*n);      // size (natom by ndim)
    float *Sk = (float *) malloc(sizeof(float)*n);      // size (natom by ndim)
    float *Yk = (float *) malloc(sizeof(float)*n);      // size (natom by ndim)
    float *Dk = (float *) malloc(sizeof(float)*n);      // size (natom by ndim)
    int *pivot = (int *) malloc(sizeof(int)*n);      // size (natom by ndim)
    // random number generator
    unsigned long int now = static_cast<unsigned long int>(time(NULL));
    MTRand *mtrand = new MTRand(now);

    // read xyz positions and 
    readInput("coordinates/Configuration.xyz", natoms, XYZ);

    // initialize vectors
    FormIdentityMatrix(n, Bjk);
    VectorFill(n, Pk, 0);
    VectorFill(n, Yk, 0);
    VectorFill(n*n, Cjk, 0);
    VectorFill(n*n, Djk, 0);

    // initialize scalars
    float ep, maxf; float c = 0; float d = 0;
    float alpha = 0.001;

    for(long iter = 0; iter < numruns; iter++) {
        long nshow = 5; // n
        printf("\n=======================");
        printf("\n__ iteration = %ld __\n", iter);   
        printf("\n=======================\n");
        if(PNT) printf("\n__XYZ_j__\n");
        if(PNT) PrintMatrix(nshow, 1, XYZ);

        /* step 0. clear all forces (n by 1) */
        ClearMatrix(n*n, 1, Cjk);
        ClearMatrix(n*n, 1, Djk);
        ClearMatrix(n*n, 1, Ejk);
        ClearMatrix(n, 1, Fj);
        ClearMatrix(n, 1, FjShift);
        
        /* step 1. compute fx, fy, fz vector (n by 1) */
        ComputeMorseForces(natoms, XYZ, Fj);

        /* step 2. solve for Pk vector */
        int err; int nrhs = 1;
        cblas_scopy(n, Fj, 1, Pk, 1);
        for(long q = 0; q < n*n; q++) mat[q] = Bjk[q]; 
        for(long q = 0; q < n; q++) vec[q] = Fj[q];
        dgesv_(&n, &nrhs, mat, &n, pivot, vec, &n, &err);
        for(long q = 0; q < n*n; q++) Bjk[q] = mat[q]; 
        for(long q = 0; q < n; q++) Pk[q] = vec[q];
        if(err > 0) printf("factorization has been completed (U is exactly singular)\n");
        if(err < 0) printf("argument had an illegal value\n");
        
        if(PNT) printf("\n__mat__\n");
        if(PNT) for(long q = 0; q < nshow; q++) printf("%ld %lf\n", q, mat[q]);
        if(PNT) printf("\n__vec__\n");
        if(PNT) for(long q = 0; q < nshow; q++) printf("%ld %lf\n", q, vec[q]);
        
        if(1) printf("\n__B_jk_start__\n");
        if(1) PrintMatrix(nshow, 1, Bjk);
        if(PNT) printf("\n__F_j__\n");
        if(PNT) PrintMatrix(nshow, 1, Fj);
        if(PNT) printf("\n__P_k__\n");
        if(PNT) PrintMatrix(nshow, 1, Pk);

        /* step 3. solve for Sk = alpha * Pk */
        VectorScalarProduct(n, alpha, Pk, Sk);
        
        if(PNT) printf("\n__S_k__\n");
        if(PNT) PrintMatrix(nshow, 1, Sk);
        
        /* step 4. solve for Xk+1 = Xk + Sk evolved position */
        VectorAddition(n, XYZ, Sk, sXYZ);
        
        if(PNT) printf("\n__XYZ_j__\n");
        if(PNT) PrintMatrix(nshow, 1, XYZ);
        if(PNT) printf("\n__sXYZ_j__\n");
        if(PNT) PrintMatrix(nshow, 1, sXYZ);

        /* step 5. solve for Yk = -FjShifted + Fj */
        ComputeMorseForces(natoms, sXYZ, FjShift);
        VectorSubtraction(n, FjShift, Fj, Yk);
        
        if(PNT) printf("\n__F_jShift__\n");
        if(PNT) PrintMatrix(nshow, 1, FjShift);
        if(1) printf("\n__Y_k__\n");
        if(1) PrintMatrix(nshow, 1, Yk);

        /* step 6. 
            solve for Bjk' = Bjk + (Yk YkT) / (YkT Sk) + (Bjk Sk SkT BjkT) / (SjT Bjk Sj)
                      Bjk' = Bjk + Cjk / c             + (Bjk Sk Bjk Sk) / (SjT Dk)
                      Bjk' = Bjk + Cjk                 + (Bjk Sk Dk) / d
                      Bjk' = Bjk + Cjk                 + (Bjk Djk) / d
                      Bjk' = Bjk + Cjk                 + Djk / d            
                      Bjk' = Bjk + Cjk                 + Djk                */
        
        /* step 6.1 solve for c (second term denominator) */
        c = VectorDotProduct(n, Yk, Sk); 
        
        if(1) printf("\n__C__ = %f\n", c);
        
        /* step 6.2 solve for Cjk (second term numerator) */
        //VectorOuterProduct(n, Yk, Yk, Cjk);
        cblas_sger(CblasRowMajor, n, n, 1, Yk, 1, Yk, 1, Cjk, n);
        VectorScalarProduct(n*n, 1/c, Cjk, Cjk);
        
        if(1) printf("\n__C_jk__\n");
        if(1) PrintMatrix(nshow, 1, Cjk);
        
        /* step 6.3 solve for Dk and then scalar d (third term denominator) */
        MatrixVectorProduct(n, Bjk, Sk, Dk);
        d = VectorDotProduct(n, Sk, Dk);
        
        if(1) printf("\n__D__ = %f\n", d);
        
        /* step 6.4 solve for Ejk (third term numerator) */
        MatrixVectorProduct(n, Bjk, Sk, Dk);
        VectorOuterProduct(n, Sk, Dk, Djk);     // here Djk is nonzero
        MatrixMatrixProduct(n, Bjk, Djk, Ejk);  // here Djk is zero
        VectorScalarProduct(n*n, 1/d, Ejk, Ejk);

        if(1) printf("\n__E_jk__\n");
        if(1) PrintMatrix(nshow, 1, Ejk);

        /* step 7. Bj+1 = Bjk + Cjk + Djk */
        float scalar = 1; long stride = 1;
        cblas_saxpy(n*n, scalar, Cjk, stride, Bjk, stride);
        cblas_saxpy(n*n, scalar, Ejk, stride, Bjk, stride);
        if(1) printf("\n__B_jk_FINAL__\n");
        if(1) PrintMatrix(nshow, 1, Bjk);

        /* step 8. potential energy summary */
        ep = ComputeTotalMorsePotential(natoms, XYZ);
        printf("\n__EP__\n");
        printf("\nEP = %f\n", ep);

        /* step 9. save shifted xyz to xyz */
        cblas_scopy(n, sXYZ, 1, XYZ, 1);

        /* step 10. maximum force of all the components */
        int idMax = GetAbsMaxElementIndex(n, 1, Fj, 1);
        maxf = Fj[idMax];

        //cout << std::scientific;
        cout.precision(2);
        cout << "\nAbs. Max force component is = " << maxf << "\n";

        Coords2XYZFile(natoms, XYZ, iter);
        if( !(n % 2) ) Results2File(iter, ep, maxf);
        if(abs(maxf) < (1e-4) ) {
            break;
        }
    }

    // free memory 
    free(vec); free(mat);
    free(XYZ); free(sXYZ);
    free(Fj); free(FjShift); 
    free(Bjk); free(Cjk); free(Djk); free(Dk); free(Ejk);
    free(Pk);  free(Sk); free(Yk); 
    free(pivot);
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
    //if(__APPLE__) {
    //    catlas_sset(n, val, v, 1);
    //} else {
        for(long i=0; i<n; i++) {
            v[i] = val;
        }
    //}
    return v;
}

// LU Decomposition
float* FormIdentityMatrix(long n, float* A) {
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

float* LUDecomposition(long n, float* bjk) {
    // recompose bjk matrix as upper triangular matrix from 
    // bjk * pk = fj equation. this function determines U matrix.
    // 'n' represents '3N' degrees of freedom.
    if(DBG) printf("LUDecomposition()\n");
    long i, j, k;
    long index1, index2;
    float c;
    for(i = 0; i < n; i++) {
        for(j = i + 1; j < n; j++) {
            c = bjk[ i*n + j ] / bjk[ i*(n+1) ];
            for(k = 0; k < n; k++) {
                //printf("in + j bjk[%ld] = %f\n", i*n + j, bjk[i*n + j]);
                //printf("i(n + 1) bjk[%ld] = %f\n", i*(n + 1), bjk[i*(n+1)]);
                //printf("jn + k, bjk[%ld] = %f\n", j*n + k, bjk[j*n + k]);
                bjk[j*n + k] = bjk[j*n + k] - bjk[i*n + k]*c;
            }
        }
    }
    return bjk;
}

float* SolvePk(long n, float* Bjk, float* Pk, float* Fj) {
    /* solve the Aij Xi = Bj linear equation. in our case we are solving for vector Pk
       Bjk Pk = -Del_j f(Xk) = Fj. Where Fj is our force vector {x1,y1,z1, x2,y2,z2,..., xn,yn,zn}

       solve for B_j^(n-1) <- superscript is prime */
    if(DBG) printf("SolvePk()\n");
    long i, j, k;
    for(i = 0; i < n; i++) {
        Pk[i] = Fj[i];
        for(j = 0; j < i; j++) {
            Pk[i] = Pk[i] - Pk[j] * Bjk[i*n + j];
            //printf("Pk[%ld] = %f Fj[%ld] = %f Bjk[%ld][%ld] = %f\n", i, Pk[i], i, Fj[i], i, j, Bjk[i*n + j]);
        }
    }

    /* now we may be able to solve for Xi or Pk vectors
        Xi = (Bi - sum_{i = j+1}^{n} Aij Xj) / Aii
        
        or

        Pk = (Fj - sum_{i = j+1}^{n} Bjk Pk) / Bjj */
    for(i = 0; i < n; i++) {
        for(j = i + 1; j < n; j++) {
            Pk[i] -= Bjk[ i*n + j ] * Pk[j];
        }
        printf("Bjk[ i*(n+1) ] = %f\n", Bjk[ i*(n+1) ]);
        Pk[i] = Pk[i] / Bjk[ i*(n+1) ];
    }
    return Pk;
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
        cout.precision(4);
        cout << setw(8) << i << "\t "<< vec[i] << "\t\n";
        //if( !((i+1)%m) ) { cout << "\n"; }
    }
}

void Coords2XYZFile(long n, float *xyz, long index) {
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

// NEW POSITIONS

float* ComputeNewPositions(long n, float* xyzvec, float* g, float* fvec) {
    // compute evolved xyzvec
    if(DBG) printf("ComputeNewPositions\n");
    for(long i = 0; i < n; i++) {
        xyzvec[i] = xyzvec[i] + g[i] * fvec[i];
    }
    return xyzvec;
}

