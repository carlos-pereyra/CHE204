#include <cstdlib>
#include <cstdio>
#include <omp.h>
//#include <mkl.h> // not installed currently
#include <vector>
#include <algorithm>
#include <iostream> // std
#include <iomanip>  // setw
#include <fstream>  // fopen
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.c
//

using namespace std;

float* MorsePotential(const long n, float *ri, float *rj, float *ep);
float* dMorsePotential(const long n, float *ri, float *rj, float *dep);
float* dMorsePotential_ForwardDiff(const long n, float *r, float *ep, float *depf);
float* dMorsePotential_BackwardDiff(const long n, float *r, float *ep, float *depf);
float* dMorsePotential_CentralDiff(const long n, float *r, float *ep, float *depf);
float* VectorFill(const long n, float *v, float val);
float* VectorSubtraction(const long n, float *v1, float *v2, float *vresult);
void write2file(const long n, float *r, float *ep, float *dep, float *dep_fd, float *dep_bd, float *dep_cd, float *err_fd, float *err_bd, float *err_cd);

//reference function to verify data
float* MorsePotential(const long n, float *r, float *ep) {
    // this function computes the morse potential.
    //
    // function takes the difference between two vectors
    // ri and rj both of length n.
    //
    // return a list containing potential energy.
    float d = 2;
    float a = 1;
    float ro = 1.2;
    for(long i=0; i<n; i++){
        ep[i] = d*pow((1 - exp(-a*(r[i]-ro))), 2);
    }
    return ep;
}

float* dMorsePotential(const long n, float *r, float *dep) {
    // compute the negative first derivative of the morse potential.
    //
    // function takes the difference between two vectors
    // ri and rj both of length n.
    //
    // return a list containing first spacial derivative of potential energy.
    float d = 2;
    float a = 1;
    float ro = 1.2;
    for(long i=0; i<n; i++){
        dep[i] = -2 * d * a * ( exp(-a*(r[i]-ro)) - exp(-2*a*(r[i]-ro)) );
    }
    return dep;
}

float* dMorsePotential_ForwardDiff(const long n, float *r, float *ep, float *depf) {
    // compute the negative first derivative of the morse potential using
    // forward difference method
    for(int i=0; i<n; i++) {
        depf[i] = -1*(ep[i+1] - ep[i]) / (r[i+1] - r[i]);
    }
    return depf;
}

float* dMorsePotential_BackwardDiff(const long n, float *r, float *ep, float *depf) {
    // compute the negative first derivative of the morse potential using
    // forward difference method
    for(int i=1; i<=n; i++) {
        depf[i] = -1*(ep[i] - ep[i-1]) / (r[i] - r[i-1]);
    }
    return depf;
}

float* dMorsePotential_CentralDiff(const long n, float *r, float *ep, float *depf) {
    // compute the negative first derivative of the morse potential using
    // forward difference method
    for(int i=1; i<n; i++) {
        depf[i] = -1*(ep[i+1] - ep[i-1]) / (r[i+1] - r[i-1]);
    }
    return depf;
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

void write2file(const long n, float *r, float *ep, float *dep, float *dep_fd, float *dep_bd, float *dep_cd, float *err_fd, float *err_bd, float *err_cd) {
    string filename = "data/output.dat";
    ofstream outfile;
    if (1) {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# r";
        outfile << left << setw(12) << "ep";
        outfile << left << setw(12) << "dep";

        outfile << left << setw(12) << "dep_fd";
        outfile << left << setw(12) << "dep_bd";
        outfile << left << setw(12) << "dep_cd";

        outfile << left << setw(12) << "err_fd";
        outfile << left << setw(12) << "err_bd";
        outfile << left << "err_cd\n";
    }
    /*
    else {
        outfile.open(filename, std::fstream::app);
    }
    */
    for(int i=0; i<=n; i++) {
        outfile << std::scientific;
        outfile.precision(4);
        outfile << left << setw(12) << r[i];
        outfile << left << setw(12) << ep[i];
        outfile << left << setw(12) << dep[i];

        outfile << left << setw(12) << dep_fd[i];
        outfile << left << setw(12) << dep_bd[i];
        outfile << left << setw(12) << dep_cd[i];

        outfile << left << setw(12) << err_fd[i];
        outfile << left << setw(12) << err_bd[i];
        outfile << left << err_cd[i] << "\n";
    }
    outfile.close();
}


int main(int argc, char** argv) {
    const long n = 101; //rows
    const long m = 1; //columns
    //float *mat = (float *) malloc((long)sizeof(float)*n*m); // 2D
    float *ri     = (float *) malloc((long)sizeof(float)*n); // 1D
    float *rj     = (float *) malloc((long)sizeof(float)*n); // 1D
    float *r      = (float *) malloc((long)sizeof(float)*n); // 1D
    float *ep     = (float *) malloc((long)sizeof(float)*n); // 1D
    float *dep    = (float *) malloc((long)sizeof(float)*n); // 1D
    float *dep_fd = (float *) malloc((long)sizeof(float)*n); // 1D
    float *dep_bd = (float *) malloc((long)sizeof(float)*n); // 1D
    float *dep_cd = (float *) malloc((long)sizeof(float)*n); // 1D
    float *err_fd = (float *) malloc((long)sizeof(float)*n); // 1D
    float *err_bd = (float *) malloc((long)sizeof(float)*n); // 1D
    float *err_cd = (float *) malloc((long)sizeof(float)*n); // 1D
    float xi = 0.5;
    float xf = 2.5;
    float dx = (xf - xi) / float(n);
  
    // perform computation
    const double t0 = omp_get_wtime();
    // fill position vectors
    for(int i=0; i<=n; i++) ri[i] = dx*i + xi;
    rj     = VectorFill(n, rj, 0);
    r      = VectorFill(n, r, 0);
    ep     = VectorFill(n, ep, 0);
    dep    = VectorFill(n, dep, 0);
    dep_fd = VectorFill(n, dep_fd, 0);
    dep_bd = VectorFill(n, dep_bd, 0);
    dep_cd = VectorFill(n, dep_cd, 0);
    err_fd = VectorFill(n, err_fd, 0);
    err_bd = VectorFill(n, err_bd, 0);
    err_cd = VectorFill(n, err_cd, 0);
    // compute values
    r      = VectorSubtraction(n, ri, rj, r);           // differential distance
    ep     = MorsePotential(n, r, ep);                  // morse potential
    dep    = dMorsePotential(n, r, dep);                // morse potential (analytic first derivative)
    dep_fd = dMorsePotential_ForwardDiff(n, r, ep, dep_fd);     // morse potential (forward diff. method)
    dep_bd = dMorsePotential_BackwardDiff(n, r, ep, dep_bd);    // morse potential (backward diff. method)
    dep_cd = dMorsePotential_CentralDiff(n, r, ep, dep_cd);     // morse potential (central diff. method)
    err_fd = VectorSubtraction(n, dep_fd, dep, err_fd); // error (dep_fd)
    err_bd = VectorSubtraction(n, dep_bd, dep, err_bd); // error (dep_bd)
    err_cd = VectorSubtraction(n, dep_cd, dep, err_cd); // error (dep_cd)
    const double t1 = omp_get_wtime();

    // Printing perf
    printf("\nrun-time: \033[42m%8.4f\033[0m (ms)\n", t1-t0);

    // write to find for plotting
    write2file(n, r,ep,dep, dep_fd,dep_bd,dep_cd, err_fd,err_bd,err_cd);

    // free mem
    free(ri);
    free(rj);
    free(r);
    free(ep);
    free(dep);
    free(dep_fd);
    free(dep_bd);
    free(dep_cd);
    free(err_fd);
    free(err_bd);
    free(err_cd);
}
