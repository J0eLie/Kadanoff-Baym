#include "propagate.h"

void Propagate::epprop(Grid &grid) {
  int nq = grid.nq;
  utk = new complex<double>[nq];
  utkc = new complex<double>[nq];
  utkconj = new complex<double>[nq];
  utkcconj = new complex<double>[nq];
  double amee = 0.067;
  double h2m2 = h2m02 / amee;
  double epsp;
  for(int iq = 1; iq < nq; iq++){
    epsp = h2m2*grid.p2[grid.idxmap[iq]];
    utk[iq] = exp(iu*epsp*dt/hbar);
    utkc[iq] = (utk[iq] - 1.0) / epsp;
    utkconj[iq] = conj(utk[iq]);
    utkcconj[iq] = conj(utkc[iq]);
  }
  // p2 = 0
  utk[0] = 1.0;
  utkc[0] = iu*dt/hbar;
  utkconj[0] = 1.0;
  utkcconj[0] = -iu*dt/hbar;

  return;
}

void Propagate::stepping(GF &gf, Collision &coll, Sigma &sig, Grid &grid, int it) {
  // First pass
  // update gf with I(T+dt) = I(T)
  int nq = gf.nq;
  int nt = gf.nt;
  double dtoverhbar = dt / hbar;
  int dim1 = nt*nq;
  // time-axis t2
  for(int it1 = 0; it1 <= it; it1++){
    for(int iq = 0; iq < nq; iq++){
      gf.gre[it1*dim1 + (it+1)*nq + iq] = utk[iq] * gf.gre[it1*dim1 + it*nq + iq] + utkc[iq]*(coll.Ile[it1*nq + iq] + coll.Ile_hf[it1*nq + iq]);
    }
  }
  // time-axis t1
  for(int it2 = 0; it2 < it; it2++){
    for(int iq = 0; iq < nq; iq++){
      gf.gre[(it+1)*dim1 + it2*nq + iq] = utkconj[iq] * gf.gre[it*dim1 + it2*nq + iq] + utkcconj[iq]*(coll.Igr[it2*nq + iq] + coll.Igr_hf[it2*nq + iq]);
    }
  }
  for(int iq = 0; iq < nq; iq++){
    gf.gre[(it+1)*dim1 + it*nq + iq] = utkconj[iq] * (gf.gre[it*dim1 + it*nq + iq] - iu) + utkcconj[iq]*(coll.Igr[it*nq + iq] + coll.Igr_hf[it*nq + iq]);
  }
  // time-axis diagonal
  for(int iq = 0; iq < nq; iq++){
    gf.gre[(it+1)*dim1 + (it+1)*nq + iq] = gf.gre[it*dim1 + it*nq + iq] + iu*dtoverhbar*(coll.Ile[it*nq + iq] - coll.Igr[it*nq + iq]);
  }

  // update self-energies
  /*sig.calSigfk(grid, gf, it+1);
  sig.calSig(grid, gf, it+1);
    
  // save collision terms at time T
  coll.savColl(it);
  // update collision terms
  coll.calColl(gf, sig, it+1);

  // Second pass
  // time-axis t2
  for(int it1 = 0; it1 <= it; it1++){
    for(int iq = 0; iq < nq; iq++){
      gf.gre[it1*dim1 + (it+1)*nq + iq] = utk[iq] * gf.gre[it1*dim1 + it*nq + iq] 
                                          + utkc[iq]*0.5*(coll.Ile[it1*nq + iq] + coll.Ilesave[it1*nq + iq]
                                          + coll.Ile_hf[it1*nq + iq] + coll.Ilesave_hf[it1*nq + iq]);
    }
  }
  // time-axis t1
  for(int it2 = 0; it2 < it; it2++){
    for(int iq = 0; iq < nq; iq++){
      gf.gre[(it+1)*dim1 + it2*nq + iq] = utkconj[iq] * gf.gre[it*dim1 + it2*nq + iq] 
                                          + utkcconj[iq]*0.5*(coll.Igr[it2*nq + iq] + coll.Igrsave[it2*nq + iq]
                                          + coll.Igr_hf[it2*nq + iq] + coll.Igrsave_hf[it2*nq + iq]);
    }
  }
  for(int iq = 0; iq < nq; iq++){
    gf.gre[(it+1)*dim1 + it*nq + iq] = utkconj[iq] * (gf.gre[it*dim1 + it*nq + iq] - iu) 
                                       + utkcconj[iq]*0.5*(coll.Igr[it*nq + iq] + coll.Igrsave[it*nq + iq]
                                       + coll.Igr_hf[it*nq + iq] + coll.Igrsave_hf[it*nq + iq]);
  }
  // time-axis diagonal
  for(int iq = 0; iq < nq; iq++){
    gf.gre[(it+1)*dim1 + (it+1)*nq + iq] = gf.gre[it*dim1 + it*nq + iq] 
                                           + iu*dtoverhbar*0.5*(coll.Ile[it*nq + iq] + coll.Ilesave[it*nq + iq] 
                                           - coll.Igr[it*nq + iq] - coll.Igrsave[it*nq + iq]);
  }*/

  return;
}

Propagate::Propagate(double dt_, Grid &grid) {
  dt = dt_;
  epprop(grid);
}

Propagate::~Propagate() {
  delete[] utk;
  delete[] utkc;
  delete[] utkconj;
  delete[] utkcconj;
}
