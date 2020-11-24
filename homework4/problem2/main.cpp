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
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.cpp
//
//      ./test
//
//      for extra juiciness use MPI - currently does not work :(
//      
//      mpicxx main.cpp -o test
//
//      mpirun -n 2 ./test
//

#define NX 80
#define NY 80
#define SCALE 6
using namespace std;

double f_lm(double x_l, double y_m);
void write2file(string filename, int i, int newlineflag, double x, double y, double z, double q);

int main(int argc, char** argv) {
    double x, y, z;
    double xmin = -SCALE*2.4; double xmax = SCALE*2.4; double ymin = -SCALE*2.4; double ymax = SCALE*2.4;
    double dx = (xmax - xmin) / NX; double dy = (ymax - ymin) / NY;
    double rx[2] = {2.4, 1.2};
    double ry[2] = {0, 2.078};
    double theta_x = 2*M_PI / NX;
    double theta_y = 2*M_PI / NY;
    double greal, gimag;
    double zreal, zimag, zabs;
    // parameters
    printf("PRECONDITIONERS\n");
    printf("===============\n");
    printf("dx = %lf\n", dx);
    printf("dy = %lf\n", dy);
    printf("M_PI = %lf\n", M_PI);

    int j, k;
    int l, m;

    // ORIGINAL FUNCTION

    for(j=0; j < NY; j++) {
        greal = 0;
        gimag = 0;
        for(k=0; k < NX; k++) {
            x = xmin + k*dx;
            y = ymin + j*dy;
            z = f_lm(x,y);
            write2file("data/function.dat", j+k, 0, x, y, z, 0);
        }
        write2file("data/function.dat", j+k, 1, x, y, z, 0);
    }

    // 2D DISCRETE FOURIER TRANSFORM
    // N_l = N_j = Nx
    // N_m = N_k = Ny
    //int myid, numprocs, n = NX;
    // we needed this to go a little faster
    //MPI_Init(&argc, &argv);
    // gets number of processors
    //MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // gets rank id number
    //MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    //MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //if(n==0) {
    //    exit(1);
    //}
    //double zr;

    #pragma omp for
    for(j=0; j < NY; j++) {
        for(k=0; k < NX; k++) {
            //zr = 0;
            zreal = 0;
            zimag = 0;
            zabs = 0;

            //====================== 2D SUM ============================
            for(l=0; l < NY; l++) {
                greal = 0;
                gimag = 0;
                for(m=0; m < NX; m++) {
                    x = xmin + m*dx;
                    y = ymin + l*dy;
                    greal +=  f_lm(x,y)*cos(j*theta_x*l + k*theta_y*m);
                    gimag += -f_lm(x,y)*sin(j*theta_x*l + k*theta_y*m);
                }
                zreal += greal;
                zimag += gimag;

            }
            //MPI_Reduce(&zreal, &zr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            // absolute value
            zabs = sqrt(pow(zreal, 2) + pow(zimag, 2));
            // normalization
            zabs /= sqrt(NX*NY);
            zreal /= sqrt(NX*NY);
            zimag /= sqrt(NX*NY);
            // print results to data file
            write2file("data/dft2d.dat", j+k, 0, j, k, zreal, zimag);
            write2file("data/dft2d_abs.dat", j+k, 0, j, k, zabs, 0);
            //============================================================
        
        }
        // print newline in data file between rows
        write2file("data/dft2d.dat", j+k, 1, j, k, zreal, zimag);
        write2file("data/dft2d_abs.dat", j+k, 1, j, k, zabs, 0);

    }
    //MPI_Finalize();

}

double f_lm(double x_l, double y_m) {
    // series function
    //
    //      f(x[l],y[m]) = sum_{n=-4}^{+4} sum_{m=-4}^{+4} exp(-(x-nrx1-mrx2)**2/sigma - (y-nry1-mry2)**2/sigma)
    //
    double fx, fy;
    double sigma = 0.1;
    double rx[2] = {2.4, 1.2};
    double ry[2] = {0, 2.078};
    double z = 0;
    int n, m;
    for(n=-4; n <= 4; n++) {
        for(m=-4; m <= 4; m++) {
            fx = exp( -1*pow((x_l-n*rx[0]-m*rx[1]),2) / sigma );
            fy = exp( -1*pow((y_m-n*ry[0]-m*ry[1]),2) / sigma );
            z += fx*fy;
        }
    }
    return z;
}

// OUTPUT

void write2file(string filename, int i, int newlineflag, double x, double y, double z, double q) {
    ofstream outfile;
    if (i==0) {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# x";
        outfile << left << setw(12) << "y";
        outfile << left << setw(12) << "z";
        outfile << left << "q\n";
    } else if(newlineflag) {
        outfile.open(filename, std::fstream::app);
        outfile << "\n";
        outfile.close();
    } else {
        outfile.open(filename, std::fstream::app);
    }
    // print values
    outfile << std::scientific;
    outfile.precision(4);
    outfile << left << setw(12) << x;
    outfile << left << setw(12) << y;
    outfile << left << setw(12) << z;
    outfile << left << q << "\n";
    outfile.close();
}
