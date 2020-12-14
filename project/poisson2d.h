#ifndef POISSON2D_H
#define POISSON2D_H

//namespace PARTICLE_DYNAMICS_NS {

class Poisson2D {
    public:
        Poisson2D(int n, int m);
        ~Poisson2D();
        void init();            // compose matrices
        
        double smooth(double*); // finite difference
        
        void writematrix2file(std::string filename, std::string mode);

        double error;

    private:
        double *uold;
        double *u;
        double *x;
        double *y;
        double *p; // charge density 
        int mem;
        int nelem, melem;
        double dx, dy, dx2, dy2;
};

//}

#endif
