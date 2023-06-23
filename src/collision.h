#ifndef collision_h
#define collision_h

#include <complex>
#include <algorithm>
#include <omp.h>
#include "gf.h"
#include "sigma.h"
#include "constants.h"
#include "mkl.h"

using namespace std;

class Collision
{
  public:
    double dt; // time step in fs
    int nt;
    int nq;
    complex<double> *Ile=nullptr; // I2^<(p,t1,T); nt x nq
    complex<double> *Igr=nullptr; // I1^>(p,T,t2); nt x nq
    complex<double> *Ilesave=nullptr; // I2^<(p,t1,T); nt x nq
    complex<double> *Igrsave=nullptr; // I1^>(p,T,t2); nt x nq
    
    complex<double> *Ile_hf=nullptr; // I2^<(p,t1,T); nt x nq
    complex<double> *Igr_hf=nullptr; // I1^>(p,T,t2); nt x nq
    complex<double> *Ilesave_hf=nullptr; // I2^<(p,t1,T); nt x nq
    complex<double> *Igrsave_hf=nullptr; // I1^>(p,T,t2); nt x nq
    
    void calColl(GF &gf, Sigma &sig, int it);
    void savColl(int it);

    Collision(double dt_, int nt_, int nq_);
    ~Collision();
};

#endif
