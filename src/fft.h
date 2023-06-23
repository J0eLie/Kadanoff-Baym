#ifndef fft_h
#define fft_h

#include <complex>
#include <omp.h>
#include "mkl.h"
#include "fftw3.h"

using namespace std;

void zfft3(complex<double> *data, const int n[]); // e^{-ikr}

void zifft3(complex<double> *data, const int n[]);// e^{ikr}

#endif
