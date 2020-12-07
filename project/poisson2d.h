#ifndef POISSON2D_H
#define POISSON2D_H

//namespace PARTICLE_DYNAMICS_NS {

class Poisson2D {
    public:
        Poisson2D(int n, int m);
        ~Poisson2D();
        // initialization
        void init();
        // computation
        void smooth(double*);
        // 
        void writematrix2file(int index, double* matrix, string filename);

    private:
        double *Uold;
        double *U;
        double *x;
        double *y;
        double *f; // charge density rhs
        int mem;
        int nelem, melem;
};

//}

#endif
