#include "sigma.h"

void Sigma::calSigfk(Grid &grid, GF &gf, int it) { // Sigma_Fk(p, T) = i int d^3p1/(2pi)^3 G^<(p1, T, T) * V_fk(p1 - p)
  int n = grid.n;
  int nfft = grid.nfft;
  int nmesh[3] = {nfft, nfft, nfft};
  complex<double> *fftbox = new complex<double>[n]();
  gf.put_into_fftbox(fftbox, grid.idxmap, it, it);
  zfft3(fftbox, nmesh);
  complex<double> *Vtmp = new complex<double>[n];
#pragma omp parallel for
  for(int i = 0; i < n; i++) Vtmp[i] = pot->Vfk[i];
  zfft3(Vtmp, nmesh);
  
  int j1, j2, j3;
#pragma omp parallel for private(j1, j2, j3)
  for(int i1 = 0; i1 < nfft; i1++){
    j1 = (nfft - i1) % nfft;
    for(int i2 = 0; i2 < nfft; i2++){
      j2 = (nfft - i2) % nfft;
      for(int i3 = 0; i3 < nfft; i3++){
        j3 = (nfft - i3) % nfft;
        fftbox[i1*nfft*nfft + i2*nfft + i3] *= Vtmp[j1*nfft*nfft + j2*nfft + j3];
      }
    }
  }
  delete[] Vtmp;

  zifft3(fftbox, nmesh);

#pragma omp parallel for
  for(int iq = 0; iq < nq; iq++) sfk[iq] = fftbox[grid.idxmap[iq]].imag();
  delete[] fftbox;
  double dp = grid.dp;
  double scal = -oneovertpic*dp*dp*dp/n;
  cblas_dscal(nq, scal, sfk, 1);

  return;
}

void Sigma::calSig(Grid &grid, GF &gf, int it) {
  int n = grid.n;
  int nfft = grid.nfft;
  int nmesh[3] = {nfft, nfft, nfft};
  complex<double> *gle = new complex<double>[n];
  complex<double> *ggr = new complex<double>[n];
  complex<double> *hp = new complex<double>[n];
  double dp = grid.dp;
  double scal = oneovertpic*dp*dp*dp/n;
  for(int itb = 0; itb <= it; itb++){
    fill_n(gle, n, complex<double>(0.0, 0.0));
    fill_n(ggr, n, complex<double>(0.0, 0.0));
    gf.put_into_fftbox(gle, grid.idxmap, itb, it);
    zfft3(gle, nmesh);
    gf.put_into_fftbox(ggr, grid.idxmap, it, itb);
    if(itb == it){
#pragma omp parallel for
      for(int iq = 0; iq < nq; iq++) ggr[grid.idxmap[iq]] -= iu;
    }
    zfft3(ggr, nmesh);
    
    int j1, j2, j3;
#pragma omp parallel for private(j1, j2, j3)
    for(int i1 = 0; i1 < nfft; i1++){
      j1 = (nfft - i1) % nfft;
      for(int i2 = 0; i2 < nfft; i2++){
        j2 = (nfft - i2) % nfft;
        for(int i3 = 0; i3 < nfft; i3++){
          j3 = (nfft - i3) % nfft;
          hp[i1*nfft*nfft + i2*nfft + i3] = gle[j1*nfft*nfft + j2*nfft + j3] * ggr[i1*nfft*nfft + i2*nfft + i3];
        }
      }
    }

    zifft3(hp, nmesh);
    cblas_zdscal(n, scal, hp, 1); // hp(p) = int d^3p1/(2*pi)^3 G^<(p1-p) * G^>(p1)
    
#pragma omp parallel for
    for(int i = 0; i < n; i++) hp[i] *= pot->V[i] * pot->V[i];
    zfft3(hp, nmesh);

#pragma omp parallel for private(j1, j2, j3)
    for(int i1 = 0; i1 < nfft; i1++){
      j1 = (nfft - i1) % nfft;
      for(int i2 = 0; i2 < nfft; i2++){
        j2 = (nfft - i2) % nfft;
        for(int i3 = 0; i3 < nfft; i3++){
          j3 = (nfft - i3) % nfft;
          gle[i1*nfft*nfft + i2*nfft + i3] *= hp[j1*nfft*nfft + j2*nfft + j3];
          ggr[i1*nfft*nfft + i2*nfft + i3] *= hp[i1*nfft*nfft + i2*nfft + i3];
        }
      }
    }

    zifft3(gle, nmesh);
#pragma omp parallel for
    for(int iq = 0; iq < nq; iq++) sle[itb*nq + iq] = gle[grid.idxmap[iq]];
    cblas_zdscal(nq, 2.0*scal, sle + itb*nq, 1);
    zifft3(ggr, nmesh);
#pragma omp parallel for
    for(int iq = 0; iq < nq; iq++) sgr[itb*nq + iq] = ggr[grid.idxmap[iq]];
    cblas_zdscal(nq, 2.0*scal, sgr + itb*nq, 1);
  }
  
  delete[] gle;
  delete[] ggr;
  delete[] hp;

  return;
}

Sigma::Sigma(int nt_, Grid &grid) {
  nt = nt_;
  nq = grid.nq;
  sfk = new double[nq];
  sle = new complex<double>[nt*nq];
  sgr = new complex<double>[nt*nq];
  pot = new Coulomb(grid);
}

Sigma::~Sigma() {
  delete[] sfk;
  delete[] sle;
  delete[] sgr;
  delete pot;
}

