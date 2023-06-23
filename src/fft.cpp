#include "fft.h"

void zfft3(complex<double> *data, const int n[]) { // e^{-ikr}
  fftw_plan p;
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());

  p = fftw_plan_dft_3d(n[0], n[1], n[2], (fftw_complex*)data, (fftw_complex*)data, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_cleanup_threads();

  return;
}

void zifft3(complex<double> *data, const int n[]) { // e^{ikr}
  fftw_plan p;
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());

  p = fftw_plan_dft_3d(n[0], n[1], n[2], (fftw_complex*)data, (fftw_complex*)data, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_cleanup_threads();

  return;
}

