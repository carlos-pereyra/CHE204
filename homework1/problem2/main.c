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
//      ./test

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
float* MorseForce(long s, long e, float *dr, float *r, float *f, float **);
float** MatrixFill(const long n, const long m, float **mat, float val);
float* VectorFill(const long n, float *v, float val);
float* VectorSubtraction(const long n, float *v1, float *v2, float *vresult);
float VectorSum(const long n, float *v);
void write2file(const long n, ManyIons* );

int main(int argc, char** argv) {
    // always works in 3D
    const long natoms = 3; const long ndim = 3;
    const long ninter = natoms*(natoms-1)/2;
    long i, j; 
    printf("there's n interactions = %ld\n\n", ninter);
    // create vectors
    ManyIons *mi = (ManyIons *) malloc((long)sizeof(ManyIons)*natoms);
    ManyIons mat;
    // initialize manyion list
    for(i=0; i<natoms; i++) {
        mi[i].r = (float*) malloc((long)sizeof(float)*ndim);
        mi[i].f = (float*) malloc((long)sizeof(float)*ndim);
        if(i==0) {
            // fill atom 1
            mi[0].r[0] = 0;
            mi[0].r[1] = 0;
            mi[0].r[2] = 0;
        }
        else if (i==1) {
            // fill atom 2
            mi[1].r[0] = 1;
            mi[1].r[1] = 1;
            mi[1].r[2] = 0;
        }
        else if (i==2) {
            // fill atom 3
            mi[2].r[0] = 0;
            mi[2].r[1] = -2;
            mi[2].r[2] = 0;
        }
    }
    // initialize force matrix
    mat.fx = (float**) malloc(sizeof(float*)*natoms);
    mat.fy = (float**) malloc(sizeof(float*)*natoms);
    mat.fz = (float**) malloc(sizeof(float*)*natoms);
    for(j = 0; j<natoms; j++) {mat.fx[j] = (float*) malloc(sizeof(float)*natoms);}
    for(j = 0; j<natoms; j++) {mat.fy[j] = (float*) malloc(sizeof(float)*natoms);}
    for(j = 0; j<natoms; j++) {mat.fz[j] = (float*) malloc(sizeof(float)*natoms);}
    MatrixFill(natoms, natoms, mat.fx, 0);
    MatrixFill(natoms, natoms, mat.fy, 0);
    MatrixFill(natoms, natoms, mat.fz, 0);
    // total system energy
    float global_ep = 0;
    // selected atom attributes
    float *ri    = NULL;
    float *rj    = NULL;
    float *fi    = NULL;
    float *fj    = NULL;
    float epi = 0;
    // exchange vectors
    float *ep    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    float *fx    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    float *fy    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    float *fz    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    float *dx    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    float *dy    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    float *dz    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    float *dr    = (float *) malloc((long)sizeof(float)*ninter); // 1D
    // perform computation
    const double t0 = omp_get_wtime();
    // compute exchange force atom by atom
    for(i=0; i<natoms; i++) {
        ri = mi[i].r;
        fi = mi[i].f;
        epi = mi[i].ep;
        // clear dx dy dz dr vector 
        VectorFill(natoms, dx, 0);
        VectorFill(natoms, dy, 0);
        VectorFill(natoms, dz, 0);
        VectorFill(natoms, dr, 0);
        VectorFill(natoms, fx, 0);
        VectorFill(natoms, fy, 0);
        VectorFill(natoms, fz, 0);
        VectorFill(natoms, ep, 0);
        //
        for(j=i+1; j<natoms; j++) {
            rj = mi[j].r;
            fj = mi[j].f;
            // fill vector of displacements between ai and aj
            dx[j] = ri[0] - rj[0];
            dy[j] = ri[1] - rj[1];
            dz[j] = ri[2] - rj[2];
            dr[j] = sqrt(pow(dx[j],2) + pow(dy[j],2) + pow(dz[j],2));
            printf("%ld -> %ld dx = %f dy = %f dz = %f dr = %f\n", i, j, dx[j], dy[j], dz[j], dr[j]);
        }
        // compute all interaction forces on ai from aj
        MorseForce(i+1, natoms, dx, dr, fx, mat.fx);
        MorseForce(i+1, natoms, dy, dr, fy, mat.fy);
        MorseForce(i+1, natoms, dz, dr, fz, mat.fz);
        // compute all interaction energy on ai from aj
        MorsePotential(i+1, natoms, dr, ep);
       
        mi[i].f[0] = VectorSum(natoms, fx);    // total x force on ai from aj
        mi[i].f[1] = VectorSum(natoms, fy);    // total y force on ai from aj
        mi[i].f[2] = VectorSum(natoms, fz);    // total z force on ai from aj
        mi[i].ep = VectorSum(natoms, ep);   // total potential on ai from aj

        global_ep += mi[i].ep;              // total ensemble energy
        printf("\n");
    }
    // print results
    cout << left << setw(12) << "i";
    cout << left << setw(12) << "x";
    cout << left << setw(12) << "y";
    cout << left << setw(12) << "z";
    cout << left << setw(12) << "fx";
    cout << left << setw(12) << "fy";
    cout << left << setw(12) << "fz";
    cout << left << "ep\n";
    for(int i=0; i<natoms; i++) {
        float x = mi[i].r[0];
        float y = mi[i].r[1];
        float z = mi[i].r[2];
        float fx = mi[i].f[0];
        float fy = mi[i].f[1];
        float fz = mi[i].f[2];
        cout << std::scientific;
        cout.precision(4);
        cout << left << setw(12) << i;
        cout << left << setw(12) << x;
        cout << left << setw(12) << y;
        cout << left << setw(12) << z;
        cout << left << setw(12) << mat.fx[i][i];
        cout << left << setw(12) << mat.fy[i][i];
        cout << left << setw(12) << mat.fz[i][i];
        cout << left << mi[i].ep << "\n";
    }
    printf("\nsummary");
    printf("\ntotal ep = %f \n", global_ep);

    // stop timer
    const double t1 = omp_get_wtime();

    // printing performance
    printf("\nrun-time: \033[42m%8.4f\033[0m (ms)\n", t1-t0);

    // free mem
    for(int i=0; i<natoms; i++) {
        free(mi[i].r);
        free(mi[i].f);
        free(mat.fx[i]);
        free(mat.fy[i]);
        free(mat.fz[i]);
    }
    free(mi);
    free(mat.fx);
    free(mat.fy);
    free(mat.fz);
    free(ep);
    free(fx);
    free(fy);
    free(fz);
    free(dx);
    free(dy);
    free(dz);
    free(dr);
  
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
            printf("%ld -> %ld dr = %f ep = %f\n", s, i, r[i], ep[i]);
    }
    return ep;
}

float* MorseForce(long s, long e, float *dr, float *r, float *f, float **matf) {
    // compute the force vector component (negative first derivative of the
    // morse potential).
    //
    // function takes the difference between two vectors dr,
    // and magnitude of displacement r.
    //
    // return a list containing force vector component (ie. fx).
    float d = 2;
    float a = 1;
    float ro = 1.2;
    for(long i=s; i<e; i++) {
        f[i] = 2 * d * a * ( exp(-2*a*(r[i]-ro)) - exp(-a*(r[i]-ro)) ) * dr[i] / r[i];
        matf[s-1][s-1] += f[i];
        matf[i][i] -= f[i];
    }
    return f;
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
