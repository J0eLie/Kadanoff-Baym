#ifndef propagate_h
#define propagate_h

#include <complex>
#include "grid.h"
#include "gf.h"
#include "sigma.h"
#include "collision.h"
#include "constants.h"

using namespace std;

class Propagate
{
  private:
    double dt;
    complex<double> *utk=nullptr; // nq
                                  // exp(i epsilon(p) dt/hbar)
    complex<double> *utkc=nullptr; // (utk - 1) / epsilon(p)
    complex<double> *utkconj=nullptr; // conjugate of utk
    complex<double> *utkcconj=nullptr; // conjugate of utkc

    void epprop(Grid &grid);

  public:
    void stepping(GF &gf, Collision &coll, Sigma &sig, Grid &grid, int it);
    Propagate(double dt_, Grid &grid);
    ~Propagate();
};
#endif
