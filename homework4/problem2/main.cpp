#include <cstdlib>    // c++
#include <cstdio>     // c++
#include <omp.h>
#include <vector>     // c++
#include <algorithm>  // c++
#include <iostream>   // c++
#include <iomanip>    // c++
#include <fstream>    // c++
#include <math.h>
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.cpp
//
//      ./test
//

#define NX 60
#define NY 60
#define SCALE 10
using namespace std;

double f_lm(double x_l, double y_m);
double series(double x, int n, int m);
void write2file(int i, int newlineflag, double x, double y, double z);

int main(int argc, char** argv) {
    double x, y;
    double xmin = -SCALE*2.4; double xmax = SCALE*2.4; double ymin = -SCALE*2.4; double ymax = SCALE*2.4;
    double dx = (xmax - xmin) / NX; double dy = (ymax - ymin) / NY;
    double rx[2] = {2.4, 1.2};
    double ry[2] = {0, 2.078};
    double theta_x = 2*M_PI / NX;
    double theta_y = 2*M_PI / NY;
    double greal, gimag;
    double hreal, himag;
    double zreal, zimag, zabs;
    // parameters
    printf("PRECONDITIONERS\n");
    printf("===============\n");
    printf("dx = %lf\n", dx);
    printf("dy = %lf\n", dy);
    printf("M_PI = %lf\n", M_PI);

    // N_l = N_j = Nx
    // N_m = N_k = Ny
    int j, k;
    for(j=0; j < NX; j++) {
        for(k=0; k < NY; k++) {
            // l, m loop
            zreal = 0;
            zimag = 0;
            zabs = 0;
            for(int l=0; l < NX; l++) {
                for(int m=0; m < NY; m++) {
                    x = xmin + l*dx;
                    y = ymin + m*dy;
                    greal = f_lm(x,y)*cos(theta_x*k*m);
                    gimag = f_lm(x,y)*sin(theta_y*k*m);
                    hreal = sin(theta_x*j*l);
                    himag = cos(theta_y*j*m);

                    zreal += greal*hreal;
                    zimag += gimag*himag;
                }
            }
            zabs = sqrt(pow(zreal, 2) + pow(zimag, 2));
            //printf("j = %d k = %d zreal = %lf\n", j, k, zreal);
            write2file(j+k, 0, j, k, zabs); // 2FT value
        }
        //printf("\n");
        write2file(j+k, 1, j, k, zabs); // next line in dat file
    }

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


void write2file(int i, int newlineflag, double x, double y, double z) {
    string filename = "data/output.dat";
    ofstream outfile;
    if (i==0) {
        outfile.open(filename, std::fstream::out);
        outfile << left << setw(12) << "# x";
        outfile << left << setw(12) << "y";
        outfile << left << "z\n";
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
    outfile << left << z << "\n";

    outfile.close();
}

