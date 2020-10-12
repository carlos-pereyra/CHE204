#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <iostream> 
#include <iomanip>
#include <fstream>  
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.c
//
//      ./test [npoints-has to be odd]
//
//  example:
//
//      ./test 111


using namespace std;

typedef struct list_struct {
    float *r; // r[0], r[1], r[2] represent x,y,z respectively
    float *f;
    float **fx;
    float **fy;
    float **fz;
    float ep;
} ManyIons;



float* MorsePotential(long s, long e, float *r, float *ep);
float** MorseForce(long s, long e, float *dr, float *r, float **fmat);
float** MatrixFill(const long n, const long m, float **mat, float val);
float* VectorFill(const long n, float *v, float val);
float* VectorSubtraction(const long n, float *v1, float *v2, float *vresult);
float VectorSum(const long n, float *v);
void PrintSummary(int i, float x, float y, float z, float fx, float fy, float fz, float ep);


int main(int argc, char** argv) {
    // integrate path of atom 2 from (1,1,0) to (2,0,0)
    const long natoms = 3; const long ndim = 3;
    float xi = 1; float yi = 1; float zi = 1;
    float xf = 2; float yf = 0; float zf = 0;
    long i, j;
    //
    //long npoints = 101;   // this has to be odd!
    long npoints = stoi(string(argv[1]));   // this has to be odd!
    long nbins = (npoints-1)/2;
    float lx = xf - xi; float ly = yf - yi; float lz = zf - zi;
    float hx = lx / npoints; float hy = ly / npoints; float hz = lz / npoints;
    // 
    float *ri    = NULL; float *rj    = NULL; float *fi    = NULL; float *fj    = NULL;
    float epi = 0; float epi_sum = 0;
    // vector: atom_
    ManyIons *mi = (ManyIons *) malloc((long)sizeof(ManyIons)*natoms);
    ManyIons mat;
    // vector: dr_, ep_(rij), f_(rij)
    float *ep    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    float *fx    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    float *fy    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    float *fz    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    float *dx    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    float *dy    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    float *dz    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    float *dr    = (float *) malloc((long)sizeof(float)*natoms); // 1D
    // vector: f_(xi)
    float *ftx   = (float *) malloc((long)sizeof(float)*npoints); // f(xi)
    float *fty   = (float *) malloc((long)sizeof(float)*npoints); // f(yi)
    float *ftz   = (float *) malloc((long)sizeof(float)*npoints); // f(zi)
    float *et   = (float *) malloc((long)sizeof(float)*npoints); // e(ri)
    
    // initialize mi
    for(i=0; i<natoms; i++) {
        mi[i].r = (float*) malloc((long)sizeof(float)*ndim);
        mi[i].f = (float*) malloc((long)sizeof(float)*ndim);
        // initial positions
        if(i==0) {
            // fill atom 1
            mi[0].r[0] = 0;
            mi[0].r[1] = 0;
            mi[0].r[2] = 0;
        }
        else if (i==1) {
            // fill atom 2
            mi[1].r[0] = xi;
            mi[1].r[1] = yi;
            mi[1].r[2] = zi;
        }
        else if (i==2) {
            // fill atom 3
            mi[2].r[0] = 0;
            mi[2].r[1] = -2;
            mi[2].r[2] = 0;
        }
    }
    // initialize mat
    mat.fx = (float**) malloc(sizeof(float*)*natoms);
    mat.fy = (float**) malloc(sizeof(float*)*natoms);
    mat.fz = (float**) malloc(sizeof(float*)*natoms);
    for(j = 0; j<natoms; j++) {mat.fx[j] = (float*) malloc(sizeof(float)*natoms);}
    for(j = 0; j<natoms; j++) {mat.fy[j] = (float*) malloc(sizeof(float)*natoms);}
    for(j = 0; j<natoms; j++) {mat.fz[j] = (float*) malloc(sizeof(float)*natoms);}
    MatrixFill(natoms, natoms, mat.fx, 0);
    MatrixFill(natoms, natoms, mat.fy, 0);
    MatrixFill(natoms, natoms, mat.fz, 0);

    // perform computation
    printf("\nn points = %ld\n", npoints);
    const double t0 = omp_get_wtime();
    for(int step=0; step<npoints; step++) {
        // clear matrix
        MatrixFill(natoms, natoms, mat.fx, 0);
        MatrixFill(natoms, natoms, mat.fy, 0);
        MatrixFill(natoms, natoms, mat.fz, 0);
        epi_sum = 0;
        for(i=0; i<natoms; i++) {
            ri = mi[i].r;
            fi = mi[i].f;

            // clear dx dy dz dr vector
            VectorFill(natoms, dx, 0);
            VectorFill(natoms, dy, 0);
            VectorFill(natoms, dz, 0);
            VectorFill(natoms, dr, 0);
            VectorFill(natoms, fx, 0);
            VectorFill(natoms, fy, 0);
            VectorFill(natoms, fz, 0);
            VectorFill(natoms, ep, 0);

            // compute exchange force on atom2 (loop through atomj)
            // i wanted to vectorize the dx, ..., dr loop but not enough time.
            for(j=i+1; j<natoms; j++) {
                // ignore atom 2
                //if(j==1) continue;
                rj = mi[j].r;
                fj = mi[j].f;

                // fill vector of displacements between ai and aj
                dx[j] = ri[0] - rj[0];
                dy[j] = ri[1] - rj[1];
                dz[j] = ri[2] - rj[2];
                dr[j] = sqrt(pow(dx[j],2) + pow(dy[j],2) + pow(dz[j],2));
                //printf("%ld->%ld dx = %f dy = %f dz = %f dr = %f\n", i, j, dx[j], dy[j], dz[j], dr[j]);
            }
       
            // compute all interaction forces on ai from aj
            MorseForce(i+1, natoms, dx, dr, mat.fx);
            MorseForce(i+1, natoms, dy, dr, mat.fy);
            MorseForce(i+1, natoms, dz, dr, mat.fz);

            // compute all interaction energy on ai from aj
            MorsePotential(i+1, natoms, dr, ep); 
        
            // add all forces and potentials on atom2
            epi = VectorSum(natoms, ep);
            epi_sum += epi;

            //
            //if(i==1) mi[1].ep = VectorSum(natoms, ep);
            printf("i = %ld x = %f y = %f z = %f fx = %f fy = %f fz = %f epi = %f\n", i, ri[0], ri[1], ri[2], mat.fx[i][i], mat.fy[i][i], mat.fz[i][i], epi_sum);
        } 

        // fill force and energy travelled array
        ftx[step] = mat.fx[1][1]; // fi[0];
        fty[step] = mat.fy[1][1]; // fi[1];
        ftz[step] = mat.fz[1][1]; // fi[2];
        et[step] = epi_sum;

        // print results
        printf("\n");
        //PrintSummary(step, mi[1].r[0], mi[1].r[1], mi[1].r[2], mat.fx[1][1], mat.fy[1][1], mat.fz[1][1], total_ep);
        
        // move second atom position
        mi[1].r[0] = mi[1].r[0] + hx; // xi
        mi[1].r[1] = mi[1].r[1] + hy; // yi
        mi[1].r[2] = mi[1].r[2] + hz; // zi
    }

    // simpsons integration for work
    float Wx = 0;
    float Wy = 0;
    float Wz = 0;
    for(int i=0; i<nbins; i++) { Wx = Wx + ( ftx[2*i] + 4*ftx[2*i+1] + ftx[2*i+2] ); }
    for(int j=0; j<nbins; j++) { Wy = Wy + ( fty[2*j] + 4*fty[2*j+1] + fty[2*j+2] ); }
    for(int k=0; k<nbins; k++) { Wz = Wz + ( ftz[2*k] + 4*ftz[2*k+1] + ftz[2*k+2] ); }
    Wx = Wx * hx/3;
    Wy = Wy * hy/3;
    Wz = Wz * hz/3;
    float W = Wx + Wy + Wz;
    
    // change in potential energy
    float dE = et[0]-et[npoints-1];
   
    // stop timer
    const double t1 = omp_get_wtime();
    // printing performance
    printf("\nargv[1] = %s", argv[1]);
    printf("\nW: %f", W);
    printf("\ndE: %f", dE);
    printf("\nrun-time: \033[42m%8.4f\033[0m (ms)\n", t1-t0);

    // free memory
    for(int i=0; i<natoms; i++) {
        free(mi[i].r);
        free(mi[i].f); 
    }
    free(mi);
    free(ep);
    free(fx);
    free(fy);
    free(fz);
    free(dx);
    free(dy);
    free(dz);
    free(dr);
    free(ftx);
    free(fty);
    free(ftz);
    free(et);
}

float* MorsePotential(long s, long e, float *r, float *ep) {
    // this function computes the morse potential.
    //
    // return a list containing potential energy.
    float d = 2;
    float a = 1;
    float ro = 1.2;
    for(long i=s; i<e; i++){
        ep[i] = d*pow((1 - exp(-a*(r[i]-ro))), 2);
    }
    return ep;
}

float** MorseForce(long s, long e, float *dr, float *r, float **fmat) {
    // compute the force vector component (negative first derivative of the
    // morse potential).
    //
    // return a list containing force vector component (ie. fx).
    float d = 2;
    float a = 1;
    float ro = 1.2;
    float f = 0;
    for(long i=s; i<e; i++) {
        f = 2 * d * a * ( exp(-2*a*(r[i]-ro)) - exp(-a*(r[i]-ro)) ) * dr[i] / r[i];
        fmat[s-1][s-1] += f;
        fmat[i][i] -= f;
    }
    return fmat;
}

float** MatrixFill(const long n, const long m, float **mat, float val) {
    // fill nxn matrix
    // return filled vector
    if(n!=m) printf("\n n is not equal to m! error..");
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            mat[i][j] = val;
        }
    }
    return mat;
}

float* VectorFill(const long n, float *v, float val) {
    // fill vector with single value
    // return filled vector
    for(int i=0; i<=n; i++) {
        v[i] = val;
    }
    return v;
}

float* VectorSubtraction(const long n, float *v1, float *v2, float *vresult) {
    // subtract vector1 by vector2.
    // return vresult
    for(int i=0; i<=n; i++) {
        vresult[i] = abs(v1[i] - v2[i]);
    }
    return vresult;
}

float VectorSum(const long n, float *v) {
    // sum all elements in v
    // return scalar sum
    float sum = 0;
    for(int i=0; i<n; i++) {
        sum += v[i];
    }
    return sum;
}

void PrintSummary(int i, float x, float y, float z, float fx, float fy, float fz, float ep) {
    // print results
    cout << left << setw(12) << "i";
    cout << left << setw(12) << "x";
    cout << left << setw(12) << "y";
    cout << left << setw(12) << "z";
    cout << left << setw(12) << "fx";
    cout << left << setw(12) << "fy";
    cout << left << setw(12) << "fz";
    cout << left << "ep\n";
    cout << std::scientific;
    cout.precision(4);
    cout << left << setw(12) << i;
    cout << left << setw(12) << x;
    cout << left << setw(12) << y;
    cout << left << setw(12) << z;
    cout << left << setw(12) << fx;
    cout << left << setw(12) << fy;
    cout << left << setw(12) << fz;
    cout << left << ep << "\n";
}
