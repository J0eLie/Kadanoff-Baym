#ifndef coulomb_h
#define coulomb_h

#include "grid.h"
#include "gf.h"
#include "constants.h"

class Coulomb
{
  public:
    double kappa=0.4/aindb; // screening parameter
    double *V=nullptr; // statically screened Coulomb potential; n
    double *Vfk=nullptr; // Unscreened Coulomb potential; n

    void potcoul(Grid &grid);
    void potcoul_fk(Grid &grid);
    void calakp(GF &gf, Grid &grid); // calculate the screening parameter from the electron distribution

    Coulomb(Grid &grid);
    ~Coulomb();
};

#endif
