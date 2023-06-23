#ifndef constants_h
#define constants_h

#include <complex>

// mathematical constants
const double pi = 3.141592653589793;
const double fpi = 12.566370614359172;
const double oneovertpic = 0.004031441804150; // 1 / (2*pi)^3
const std::complex<double> iu(0.0, 1.0);

// physical constants
const double e2 = 14.3997; // eV*Ang
const double hbar = 0.658212; // eV*fs
const double h2m02 = 3.80998; // eV*Ang*Ang

const double epsilon_b = 12.9980;
const double eryd = 0.42e-2; // one exciton Rydberg (Ry) in eV
const double aindb = 132; // Bohr radius (aB) of exciton in Ang

#endif
