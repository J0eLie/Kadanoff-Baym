#ifndef sigma_h
#define sigma_h

#include <complex>
#include <algorithm> // fill_n
#include "grid.h"
#include "gf.h"
#include "coulomb.h"
#include "constants.h"
#include "fft.h"
#include "mkl.h"

using namespace std;

class Sigma
{
  public:
    Coulomb *pot=nullptr; // Coulomb potential

    int nt;
    int nq;
    double *sfk=nullptr; // Sigma^Fock(p,T); nq
    complex<double> *sle=nullptr; // Sigma^<(p,t1,T); nt x nq
    complex<double> *sgr=nullptr; // Sigma^>(p,T,t2); nt x nq

    void calSigfk(Grid &grid, GF &gf, int it);
    void calSig(Grid &grid, GF &gf, int it);

    Sigma(int nt_, Grid &grid);
    ~Sigma();
};

#endif
