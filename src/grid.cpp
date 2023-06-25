#include "grid.h"

Grid::Grid(int nfft_, int idim_, double dp_) {
  nfft = nfft_;
  n = nfft*nfft*nfft;
  idim = idim_;
  dp = dp_;
  int tmp[3];
  int *idxmap_t = new int[n];
  double *mesh_t = new double[n*3];
  int idx = 0;
  int iq = 0;
  p2 = new double[n];
  for(int i = 0; i < nfft; i++){
    tmp[0] = i > nfft/2 ? i - nfft : i;
    for(int j = 0; j < nfft; j++){
      tmp[1] = j > nfft/2 ? j - nfft : j;
      for(int k = 0; k < nfft; k++){
        tmp[2] = k > nfft/2 ? k - nfft : k;
        p2[idx] = tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2];
        if(p2[idx] <= idim*idim){
          for(int ii = 0; ii < 3; ii++) mesh_t[iq*3 + ii] = tmp[ii]*dp;
          idxmap_t[iq] = idx;
          iq += 1;
        }
        p2[idx] *= dp*dp;
        idx += 1;
      }
    }
  }
  
  nq = iq;
  idxmap = new int[nq];
#pragma omp parallel for
  for(int iq = 0; iq < nq; iq++) idxmap[iq] = idxmap_t[iq];
  delete[] idxmap_t;
  mesh = new double[nq*3];
#pragma omp parallel for
  for(int iq = 0; iq < nq; iq++){
    for(int i = 0; i < 3; i++){
      mesh[iq*3 + i] = mesh_t[iq*3 + i];
    }
  }
  delete[] mesh_t;

}

Grid::~Grid() {
  delete[] idxmap;
  delete[] mesh;
  delete[] p2;
}
