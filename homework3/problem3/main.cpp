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
#include <algorithm>    // std::sort
#include <vector>       // std::vector

#define DBG 0
#define PNT 0

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

using namespace std;

// INITIALIZATION
double* ReadInput(string filename, long nlines, double* mat);

// IO
void PrintMatrix(string desc, int m, int n, double* a, int lda );
void WriteResults(long n, double* ei, double* vi, double* ei2, double* vi2 );

// HESSIAN
double* ComputeHessian(long n, double* x, double* hij);

#define NATOMS 32    // natoms matches the number of atoms in .xyz file
#define NDIM 3

#define N NATOMS*NDIM
#define LDA N

#define NRUNS 1100
#define NSHOW 5

#define ALPHA 0.0001
#define MASS 23
#define STRIDE 1

int main(int argc, char** argv) {
    /*
     * VECTOR & MATRICES
     */
    int mem =0;
    double *Xi    = (double *) malloc(sizeof(double)*N); mem += sizeof(double)*N;     // N x 1
    double *Xi_2  = (double *) malloc(sizeof(double)*N); mem += sizeof(double)*N;     // N x 1
    
    double *Kij   = (double *) malloc(sizeof(double)*N*N); mem += sizeof(double)*N*N; // N x N
    double *Kij_2 = (double *) malloc(sizeof(double)*N*N); mem += sizeof(double)*N*N; // N x N
    
    double* Ei    = (double *) malloc(sizeof(double)*N); mem += sizeof(double)*N;     // N x 1
    double* Vi    = (double *) malloc(sizeof(double)*3*N); mem += sizeof(double)*3*N; // N x 3
    
    double* Ei_2  = (double *) malloc(sizeof(double)*N); mem += sizeof(double)*N;     // N x 1
    double* Vi_2  = (double *) malloc(sizeof(double)*3*N); mem += sizeof(double)*3*N; // N x 3

    char JOBZ='V'; char UPLO='U'; int N_INT = N; int LDA_INT = LDA; int LWORK = 3*N - 1; int INFO;
    
    /* 1.1 READ COORDINATES (MY OPTIMIZED GEOMETRY) */
    ReadInput("coordinates/Final_Coords.xyz", NATOMS, Xi);
    //if(1) PrintMatrix("Xi", 1, 10, Xi, 1);

    /* 1.2 COMPUTE HESSIAN Kij = d''V / dxi dxj */
    ComputeHessian(N, Xi, Kij);
    PrintMatrix("Hessian", N, N, Kij, N);
    
    /* 1.3 SCALE Kij BY MASS */
    cblas_dscal(N*N, 1/ (double) MASS, Kij, STRIDE); // Kij = Kij / M

    /* 1.4 FIND EIGENVALUES OF Kij */
    dsyev_(&JOBZ, &UPLO, &N_INT, Kij, &LDA_INT, Ei, Vi, &LWORK, &INFO);
   

    /* 2.1 READ COORDINATES (TRUE OPTIMIZED GEOMETRY) */
    //ReadInput("coordinates/Optimized_Argon_Cluster_Morse.xyz", NATOMS, Xi_2);
    
    /* 2.2 COMPUTE HESSIAN Kij = d''V / dxi dxj */
    //ComputeHessian(N, Xi_2, Kij_2);
    
    /* 2.3 SCALE Kij BY MASS */
    //cblas_dscal(N*N, 1/ (double) MASS, Kij_2, STRIDE); // Kij = Kij / M
    
    /* 2.4 FIND EIGENVALUES OF Kij */
    //dsyev_(&JOBZ, &UPLO, &N_INT, Kij_2, &LDA_INT, Ei_2, Vi_2, &LWORK, &INFO);
   

    /*
     * PRINT / WRITE RESULTS
     */
    WriteResults(N, Ei, Vi, Ei_2, Vi_2);
    PrintMatrix("Eigenvalues", 1, N, Ei, 1);
    //PrintMatrix("Eigenvalues Optimized", 1, N, Ei_2, 1);

    /*
     * FREE MEMORY
     */
    printf("\nfree %lf (kB)\n", mem / 1000.);

    free(Xi);   free(Kij);   free(Ei);   free(Vi);
    free(Xi_2); free(Kij_2); free(Ei_2); free(Vi_2);
}

//
// INITIALIZATION
//

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
    while(std::getline(infile, line) && counter<nlines) {
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

//
// IO
//

void PrintMatrix(string desc, int m, int n, double* a, int lda ) {
        int i, j, c = 0;
        std::cout << "\n" << desc << "\n";
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) {
                    std::cout << setprecision(2);
                    std::cout << setw(12) << a[i*lda+j]; c++;
                    //printf( "%d -> %6.6f\n", c,a[i*lda+j] ); c++;
                }
                std::cout << "\n";
        }
}

void WriteResults(long n, double* ei, double* vi, double* ei2, double* vi2) {
    // write first set of eigen modes and second set of eigen values
    // first set corresponds to coordinates I generated.
    // second set corresponds to coordinates generated by professor.
    string filename = "data/output.dat";
    ofstream file;
    
    // open new file
    file.open(filename, std::fstream::out);
    file << left << setw(20) << "# eig-val-1";
    file << left << setw(20) << "||eig-vec-1||";
    file << left << setw(20) << "eig-val-2";
    file << "||eig-vec-2||" << "\n";

    // write results to file
    for(int i = 0; i < n; i++) {
        //file << std::scientific;
        file << left << setw(20) << ei[i];
        file << left << setw(20) << sqrt(pow(vi[i*3 + 0], 2) + pow(vi[i*3 + 1], 2) + pow(vi[i*3 + 2], 2));
        file << left << setw(20) << ei2[i];
        file << sqrt(pow(vi2[i*3 + 0], 2) + pow(vi2[i*3 + 1], 2) + pow(vi2[i*3 + 2], 2)) << "\n";
    }
    file.close();
}

//
// HESSIAN
//
double MorseForce(double xij, double rij) {
    //
    if(DBG) printf("\nMorseForce()\n");
    double d = 1; double a = 1; double ro = 3.5;
    double f = 2*d*a * ( exp(-2*a*(rij-ro)) - exp(-a*(rij-ro)) ) * xij / rij;
    return f;
}

double* ComputeMorseForcesOnEachParticle(int natom, double *x, double *f) {
    if(1) printf("\nComputeMorseForcesOnAParticle()\n");
    double d = 1; double a = 1; double ro = 3.5;
    double force; double force_sum;
    double xij, yij, zij, rij, sep;
    long i, j, k;
    // upper diagonal forces in matrix
    for(i = 0; i < natom; i++) {
        force_sum = 0;
        for(j = i + 1; j < natom; j++) {
            xij = x[i*3 + 0] - x[j*3 + 0];
            yij = x[i*3 + 1] - x[j*3 + 1];
            zij = x[i*3 + 2] - x[j*3 + 2];
            rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));
            for(k = 0; k < 3; k++) {
                sep = x[i*3 + k] - x[j*3 + k];
                force = 2*d*a * ( exp(-2*a*(rij-ro)) - exp(-a*(rij-ro)) ) * sep / rij;
                // 
                //printf("particle %ld interacts w %ld in %ld dim, force = %f\n", i, j, k, force);
                f[3*i + k] += force;
                f[3*j + k] -= force;
            }
        }
    }
    return f;
}

double* ComputeHessian(long n, double* x, double* hij) {
    /* COMPUTE N x N REAL HESSIAN
     *
     * finite central difference method of second order derivative requires differential distance
     *
     * CENTRAL SECOND DERIVATIVE
     *
     *          Kij = (-Fi(xi+h) - -Fi(xi-h)) / dx
     * 
     */
    long i, j, k, l; long count = 0;
    double xi, xj;
    double dx = ALPHA;
    double xij, yij, zij;
    double rij_plus;
    double rij_minus;
    /*
     * F(x + dx) = f1 (force)
     * F(x - dx) = f2 (force)
     */
    long natom = n / 3;
    double *xplus  = (double *) malloc(sizeof(double)*natom*3); // N x 1
    double *xminus = (double *) malloc(sizeof(double)*natom*3); // N x 1
    double *f      = (double *) malloc(sizeof(double)*natom*3); // N x 1
    double *fplus  = (double *) malloc(sizeof(double)*natom*3); // N x 1
    double *fminus = (double *) malloc(sizeof(double)*natom*3); // N x 1
    // fill xplus and xminus vectors
    for(i = 0; i < n; i++) {
        xplus[i] = x[i];
        xminus[i] = x[i];
    }
    catlas_dset(n, 0, f, STRIDE);
    ComputeMorseForcesOnEachParticle(natom, x, f);
    
    // fill hessian
    for(i = 0; i < n; i++) {
        // wiggle i'th atom
        xplus[i]  = x[i] + dx;
        xminus[i] = x[i] - dx;
        catlas_dset(n, 0, fplus, STRIDE);
        catlas_dset(n, 0, fminus, STRIDE);
        ComputeMorseForcesOnEachParticle(natom, xplus, fplus);
        ComputeMorseForcesOnEachParticle(natom, xminus, fminus);
        for(j = i; j < n; j++) {
            hij[i*n + j] = ( -fplus[j] + fminus[j] ) / (2*dx);
            hij[j*n + i] = ( -fplus[j] + fminus[j] ) / (2*dx);
            printf("\n * A[%ld][%ld]\n", i, j );
            printf(" * Total count = %ld\n", count);
            printf(" * fplus = %lf fminus = %lf", fplus[j], fminus[j] );
            printf(" * -fminus + fplus = %lf\n", -fminus[j] + fplus[j] );
            printf(" * hij = %f\n", hij[count] );
            // reset j'th atom wiggle
            xplus[j] = x[j];
            xminus[j] = x[j];
            count++;
            //}
        }
        //}
        // reset i'th atom
        xplus[i]  = x[i];
        xminus[i] = x[i];
    }
    free(xplus); free(xminus); free(f); free(fplus); free(fminus);
    return hij;
}
