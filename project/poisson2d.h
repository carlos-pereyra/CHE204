#ifndef POISSON2D_H
#define POISSON2D_H

//namespace PARTICLE_DYNAMICS_NS {

class Poisson2D {
    public:
        Poisson2D(int n, int m, double*, double*);
        ~Poisson2D();
        void init(double*, double*);                            // compose matrices
        
        void smooth(int, int, double*, double*);                // finite difference operator
        void defect(int, int, double*, double*, double*);       // defect operator
        void restriction(int, int, double*, double*);           // restriction operator
        void prolongation(int, int, double*, double*, double*); // prolongation operator

        void writematrix2file(std::string filename, std::string mode);
        void copy(int, double*, double*);

        int index;
        double error;
        double errortmp;

    private:
        double *utmp;
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
