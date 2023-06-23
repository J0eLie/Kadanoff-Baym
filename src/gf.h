#ifndef gf_h
#define gf_h

#include <complex>
#include <cmath>
#include "grid.h"
#include "constants.h"

using namespace std;

class GF
{
  // Lesser Green function at time t1 and t2
  // G^<(p, t1, t2) = i <a_p^dagger(t2)a_p(t1)>
  // Greater Green function at time t1 and t2
  // G^>(p, t1, t2) = -i <a_p(t1)a_p^dagger(t2)>
  public:
    int nq; // number of momentum grid
    int nt; // number of time grid
    complex<double> *gre=nullptr; // nt x nt x nq, row-major
                                  // gre[it1, it2, iq] = G^< , it1 <= it2
                                  //                   = G^> , it1 > it2
    void init_elec(Grid &grid);
  
    complex<double> gle(int iq, int it1, int it2);
    complex<double> ggr(int iq, int it1, int it2);

    void put_into_fftbox(complex<double> *fftbox, int *idxmap, int it1, int it2);
    
    GF(Grid &grid, int nt_);
    ~GF();
};

#endif
