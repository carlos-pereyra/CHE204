#include <cstdlib>    // c++
#include <cstdio>     // c++
#include <omp.h>
#include <vector>     // c++
#include <algorithm>  // c++
#include <iostream>   // c++
#include <iomanip>    // c++
#include <fstream>    // c++
#include <math.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif
// compile this code with,
//
//      g++ -Xpreprocessor -fopenmp -lomp -o test main.c
//
//      ./test
//

#define ranfloat(w) ( w * ((rand() / (double) RAND_MAX)  - 0.5) )
#define DBG 0
using namespace std;

float* readInput(string filename, long nlines, float *mat);
float* setSMatrix(long n, long m, float sxy, float syz, float szy);
float* setXLatticePosition(long n, float *x, float dx);
float* setYLatticePosition(long n, float *y, float dy);
float* setZLatticePosition(long n, float *z, float dz);
void Coords2XYZFile(const long n, float result, float mcresult);

int main(int argc, char** argv) {
    long natoms = 30; long ndim = 3;
    // initialize
    float *xyzmat = new float(natoms*ndim);
    float *xvec = new float(natoms);
    float *yvec = new float(natoms);
    float *zvec = new float(natoms);
    // begin
    readInput("ConfigurationMod.xyz", natoms, xyzmat);
    // matrix multiply

    // clean
    free(xyzmat); free(xvec); free(yvec); free(zvec);
}

float* readInput(string filename, long nlines, float *mat) {
    // read in positions from xyz file format
    ifstream infile(filename);
    string line, elem;
    float x, y, z;
    long counter;
    if(!infile.is_open()) abort();
    counter=0;
    //while(getline(infile, line, '\n')) {
    while(counter < nlines) {
        infile >> elem >> x >> y >> z;
        cout << counter << " elem = " << elem << " x = " << x << "\n";
        counter++;
        //printf("elem = %s x = %f y = %f z = %f\n", elem.c_str(), x, y, z);
    }
    infile.close();
    return mat;
}

float* setXLatticePosition(long n, float *x, float dx) {
    // fcc lattice construction
    long i, ai; long counter = 0; long unitcells = ceil(n/4);

    for(i = 0; i < unitcells; i++) {
        for(ai=0; ai<4; ai++) {
            // atom arrangement
            if(ai == 0) x[counter] = dx * i;
            if(ai == 1) x[counter] = dx * (i+0.5);
            if(ai == 2) x[counter] = dx * (i+0.5);
            if(ai == 3) x[counter] = dx * i;
            counter++;
            if(counter>n) break;
        }
        if(counter>n) break;
    }
    return x;
}

float* setYLatticePosition(long n, float *y, float dy) {
    // fcc lattice construction
    long j, aj; long counter = 0; long unitcells = ceil(n/4);

    for(j = 0; j < unitcells; j++) {
        for(aj=0; aj<4; aj++) {
            // atom arrangement
            if(aj == 0) y[counter] = dy * j;
            if(aj == 1) y[counter] = dy * (j+0.5);
            if(aj == 2) y[counter] = dy * j;
            if(aj == 3) y[counter] = dy * (j+0.5);
            counter++;
            if(counter>n) break;
        }
        if(counter>n) break;
    }
    return y;
}

float* setZLatticePosition(long n, float *z, float dz) {
    // fcc lattice construction
    long k, ak; long counter = 0; long unitcells = ceil(n/4);

    for(k = 0; k < unitcells; k++) {
        for(ak=0; ak<4; ak++) {
            // atom arrangement
            if(ak == 0) z[counter] = dz * k;
            if(ak == 1) z[counter] = dz * k;
            if(ak == 2) z[counter] = dz * (k+0.5);
            if(ak == 3) z[counter] = dz * (k+0.5);
            counter++;
            if(counter>n) break;
        }
        if(counter>n) break;
    }
    return z;
}

int Coords2XYZFile(long n, float *x, float *y, float *z, int index) {
    // write coordinates to xyz file
    string filename_xyz = "data/output.xyz";
    ofstream myfile;
    long i; 
    if (index == 0) myfile.open(filename_xyz, std::fstream::out);
    else myfile.open(filename_xyz, std::fstream::app);

    myfile << n << "\n";
    myfile << "(n=" << index << ")!\n";
    for (i=0; i<n; i++) {
        // write atom info to file.
        myfile << "Na" << "\t";
        myfile << x[i] << "\t";
        myfile << y[i] << "\t";
        myfile << z[i] << "\n";
    }
    myfile.close();
    return 1;
}
