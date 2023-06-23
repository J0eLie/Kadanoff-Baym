#include "collision.h"

void Collision::calColl(GF &gf, Sigma &sig, int it) { // calculate the collision terms
                                                      // The trapezoidal rule is used to perform the time integrals
  // Ile
  // Hartree-Fock part G^<(p, t1, T)*Sigma^HF(p, T)
  for(int it1 = 0; it1 <= it; it1++){
    for(int iq = 0; iq < nq; iq++)
      Ile_hf[it1*nq + iq] = gf.gle(iq, it1, it) * sig.sfk[iq];
  }
  
  double dtoverhbar = dt/hbar;
  
  fill_n(Ile, (it+1)*nq, complex<double>(0.0, 0.0));
  // int_t0^t1 dt/hbar [G^>(p, t1, t) - G^<(p, t1, t)]*Sigma^<(p, t, T)
  for(int it1 = 1; it1 <= it; it1++){
    for(int iq = 0; iq < nq; iq++){
      Ile[it1*nq + iq] += 0.5*dtoverhbar*(gf.ggr(iq, it1, 0) - gf.gle(iq, it1, 0))*sig.sle[iq];
      Ile[it1*nq + iq] += 0.5*dtoverhbar*(gf.ggr(iq, it1, it1) - gf.gle(iq, it1, it1))*sig.sle[it1*nq + iq];
      for(int itb = 1; itb < it1; itb++){
        Ile[it1*nq + iq] += dtoverhbar*(gf.ggr(iq, it1, itb) - gf.gle(iq, it1, itb))*sig.sle[itb*nq + iq];
      }
    }
  }

  // - int_t0^T dt/hbar G^<(p, t1, t)*(Sigma^>(p, t, T) - Sigma^<(p, t, T))
  if(it > 0){
    for(int it1 = 0; it1 <= it; it1++){
      for(int iq = 0; iq < nq; iq++){
        Ile[it1*nq + iq] += 0.5*dtoverhbar*gf.gle(iq, it1, 0)*(conj(sig.sgr[iq]) + sig.sle[iq]);
        Ile[it1*nq + iq] += 0.5*dtoverhbar*gf.gle(iq, it1, it)*(conj(sig.sgr[it*nq + iq]) + sig.sle[it*nq + iq]);
        for(int itb = 1; itb < it; itb++){
          Ile[it1*nq + iq] += dtoverhbar*gf.gle(iq, it1, itb)*(conj(sig.sgr[itb*nq + iq]) + sig.sle[itb*nq + iq]);
        }
      }
    }
  }

  // Igr
  // Sigma^HF(p, T) * G^>(p, T, t2)
  for(int it2 = 0; it2 <= it; it2++){
    for(int iq = 0; iq < nq; iq++)
      Igr_hf[it2*nq + iq] = gf.ggr(iq, it, it2) * sig.sfk[iq];
  }

  fill_n(Igr, (it+1)*nq, complex<double>(0.0, 0.0));
  // int_t0^T dt/hbar (Sigma^>(p, T, t) - Sigma^<(p, T, t))*G^>(p, t, t2)
  if(it > 0){
    for(int it2 = 0; it2 <= it; it2++){
      for(int iq = 0; iq < nq; iq++){
        Igr[it2*nq + iq] += 0.5*dtoverhbar*(sig.sgr[iq] + conj(sig.sle[iq]))*gf.ggr(iq, 0, it2);
        Igr[it2*nq + iq] += 0.5*dtoverhbar*(sig.sgr[it*nq+iq] + conj(sig.sle[it*nq+iq]))*gf.ggr(iq, it, it2);
        for(int itb = 1; itb < it; itb++){
          Igr[it2*nq + iq] += dtoverhbar*(sig.sgr[itb*nq+iq] + conj(sig.sle[itb*nq+iq]))*gf.ggr(iq, itb, it2);
        }
      }
    }
  }

  // - int_t0^t2 dt/hbar Sigma^>(p, T, t)*(G^>(p, t, t2) - G^<(p, t, t2))
  for(int it2 = 1; it2 <= it; it2++){
    for(int iq = 0; iq < nq; iq++){
      Igr[it2*nq + iq] -= 0.5*dtoverhbar*sig.sgr[iq]*(gf.ggr(iq, 0, it2) - gf.gle(iq, 0, it2));
      Igr[it2*nq + iq] -= 0.5*dtoverhbar*sig.sgr[it2*nq + iq]*(gf.ggr(iq, it2, it2) - gf.gle(iq, it2, it2));
      for(int itb = 1; itb < it2; itb++){
        Igr[it2*nq + iq] -= dtoverhbar*sig.sgr[itb*nq + iq]*(gf.ggr(iq, itb, it2) - gf.gle(iq, itb, it2));
      }
    }
  }

  return;
}

void Collision::savColl(int it) {
  cblas_zcopy((it+1)*nq, Ile, 1, Ilesave, 1);
  cblas_zcopy((it+1)*nq, Igr, 1, Igrsave, 1);
  
  cblas_zcopy((it+1)*nq, Ile_hf, 1, Ilesave_hf, 1);
  cblas_zcopy((it+1)*nq, Igr_hf, 1, Igrsave_hf, 1);

  return;
}

Collision::Collision(double dt_, int nt_, int nq_) {
  dt = dt_;
  nt = nt_;
  nq = nq_;
  Ile = new complex<double>[nt*nq];
  Igr = new complex<double>[nt*nq];
  Ilesave = new complex<double>[nt*nq];
  Igrsave = new complex<double>[nt*nq];
  Ile_hf = new complex<double>[nt*nq];
  Igr_hf = new complex<double>[nt*nq];
  Ilesave_hf = new complex<double>[nt*nq];
  Igrsave_hf = new complex<double>[nt*nq];
} 

Collision::~Collision() {
  delete[] Ile;
  delete[] Igr;
  delete[] Ilesave;
  delete[] Igrsave;
  delete[] Ile_hf;
  delete[] Igr_hf;
  delete[] Ilesave_hf;
  delete[] Igrsave_hf;
}
