#include "gf.h"

complex<double> GF::gle(int iq, int it1, int it2) const {
  if(it1 <= it2) return gre[it1*nt*nq + it2*nq + iq];
  return -conj(gre[it2*nt*nq + it1*nq + iq]);
}

complex<double> GF::ggr(int iq, int it1, int it2) const {
  if(it1 > it2) return gre[it1*nt*nq + it2*nq + iq];
  if(it1 == it2) return -iu + gre[it1*nt*nq + it2*nq + iq];
  return -conj(gre[it2*nt*nq + it1*nq + iq]);
}

void GF::put_into_fftbox(complex<double> *fftbox, int *idxmap, int it1, int it2) { // fftbox : nfft x nfft x nfft, row-major
                                                                                   // should be filled with zeros before
#pragma omp parallel for
  for(int iq = 0; iq < nq; iq++){
    fftbox[idxmap[iq]] = gre[it1*nt*nq + it2*nq + iq];
  }

  return;
}

void GF::init_elec(Grid &grid) { // G^<(p, 0, 0) = i*f(p,0)
  double ra = 1.2;
  double rhb = 0.05;
  double rgam = 0.0194;
  double rgam1 = 6.85;
  double rgam2 = 2.1;

  gre[0] = iu*2.0/3.0*ra*exp(-0.5*rhb*rhb/rgam/rgam);
  double cost2;
  double p2;
  double fh;
  double fl;
  double chie;
  double chih;
  double chil;
  for(int iq = 1; iq < nq; iq++){
    p2 = grid.p2[grid.idxmap[iq]];
    cost2 = grid.mesh[iq*3 + 2] / sqrt(p2);
    cost2 = cost2 * cost2;
    chie = h2m02/amee * p2;
    chih = h2m02*(rgam1 - 2.0*rgam2)*p2;
    chil = h2m02*(rgam1 + 2.0*rgam2)*p2;
    fh = exp(-0.5 * (rhb - chie - chih)*(rhb - chie - chih) / rgam / rgam);
    fl = exp(-0.5 * (rhb - chie - chil)*(rhb - chie - chil) / rgam / rgam);
    gre[iq] = iu/6.0*ra*(3.0*fh*(1.0 - cost2) + fl*(1 + 3.0*cost2));
  }

  return;
}

GF::GF(Grid &grid, int nt_) {
  nt = nt_;
  nq = grid.nq;
  gre = new complex<double>[nt*nt*nq];
  init_elec(grid);
}

GF::~GF() {
  delete[] gre;
}
