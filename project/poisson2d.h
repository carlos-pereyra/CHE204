#ifndef POISSON2D_H
#define POISSON2D_H

class Poisson2D {
    public:
        Poisson2D(int n, int m, double*, double*);
        ~Poisson2D();
        void init(double*, double*);            // compose matrices
        
        void smooth(int l, int nu, double*, double*);          // finite difference operator
        void defect(int l, int nu, double*, double*, double*); // defect operator
        void restriction(double*, double*);     // restriction operator
        void prolongation(double*, double*);    // prolongation operator

        void writematrix2file(std::string filename, std::string mode);
        void write2file(std::string filename, double x, double y, double z);
        
        double error;

    private:
        double *u;
        double *x;
        double *y;
        double *f; // charge density 
        int mem;
        int nelem, melem;
        double dx, dy, dx2, dy2;
};

#endif
