#include <cstdlib>    // c++
#include <cstdio>     // c++
#include <omp.h>
#include <vector>     // c++
#include <algorithm>  // c++
#include <iostream>   // c++
#include <iomanip>    // c++
#include <fstream>    // c++
#include <math.h>
#include <sstream>
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

double** inputsignal(string filename, double **addr, double *xj, double *yj, int *n);
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

    int n = 1;
    double **address = (double**) malloc(sizeof(double*)*2);
    double *xj = (double*) malloc(sizeof(double)*n);
    double *yj = (double*) malloc(sizeof(double)*n);
    
    // input signal
    address = inputsignal("noisyspectrum.dat", address, xj, yj, &n);
    xj = address[0];
    yj = address[1];

    // parameters
    printf("PARAMETERS\n");
    printf("===============\n");
    printf("dx = %lf\n", dx);
    printf("dy = %lf\n", dy);
    printf("M_PI = %lf\n", M_PI);
    printf("n = %d\n", n);

    for(int k=0; k<n; k++) {
        //printf("xj[%d] = %lf yj[%d] = %lf\n", k, xj[k], k, yj[k]);
    }
    // N_l = N_j = Nx
    // N_m = N_k = Ny
    /*
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
    */

    // FREE MEMORY
    free(address[0]);
    free(address[1]);
    free(address);
    //free(xj);
    //free(yj);
}

// INPUT

double** inputsignal(string filename, double **addr, double *xj, double *yj, int *n) {
    // read in positions from .xyz file
    ifstream infile(filename);
    string line;
    double x, y;
    int index = 0;

    // IO
    if(!infile.is_open()) exit(1);
    //std::getline(infile, line); std::istringstream iss(line);
    //std::cout << "natoms = " << natoms << "\n";
    //std::getline(infile, line);  // line 2

    // READ LINES
    while(std::getline(infile, line)) {
        // read x y z coords
        std::istringstream iss(line);

        if (iss >> x >> y) // cout << "line = "<< index << " x = " << x << " y = " << y << "\n";

        xj[index] = double(x);
        yj[index] = double(y);
        
        index++;
        double* tmp_xj = (double* ) realloc(xj, (index+1)*sizeof(double));
        double* tmp_yj = (double* ) realloc(yj, (index+1)*sizeof(double));
        if(tmp_xj != NULL) xj = tmp_xj;
        if(tmp_yj != NULL) yj = tmp_yj;
    }
    infile.close();
    
    addr[0] = xj;
    addr[1] = yj;
    n[0] = index;

    return addr;
}

// OUTPUT

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

