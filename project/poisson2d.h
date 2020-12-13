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
        void writematrix2file(std::string filename, std::string mode);

    private:
        double *uold;
        double *u;
        double *x;
        double *y;
        double *p; // charge density 
        int mem;
        int nelem, melem;
};

//}

#endif
