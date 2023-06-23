#ifndef output_h
#define output_h

#include <iostream>
#include <iomanip>
#include <ctime>
#include <omp.h>
using namespace std;

#include "grid.h"
#include "gf.h"
#include "sigma.h"
#include "collision.h"
#include "constants.h"

class Output
{
  int outfreq;
  double dt;

  double t;
  double density;
  double Qm;
  double TEn;
  double KEn;
  double PEn_mf;
  double PEn_corr;

  public: 
    double write_header();
    void write_line(int it, const Grid &grid,  const GF &green, const Sigma &sig, const Collision &coll);
    void write_distrib(int it);
    void write_ender(double start_time);

    Output(double dt_, int pfreq);
    ~Output();
};

#endif
