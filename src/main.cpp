#include "grid.h"
#include "gf.h"
#include "sigma.h"
#include "collision.h"
#include "propagate.h"

#include <complex>
#include <iostream>
#include <string>
using namespace std;

int main(int argc, char *argv[])
{
  // ======================================================================== //
  //                                PARAMETERS                                //
  // ======================================================================== //
  if(argc != 2){
    printf("Invalid number of command line arguments\n");
    return 0;
  }
  
  // momentum grid
  int nfft = 60;
  int idim = 28;
  double dp = 0.0022;
  
  // time grid
  int ntimes = atoi(argv[1]);
  double dt = 4.0;

  // ======================================================================== //
  //                           INITIALIZE VARIABLES                           //
  // ======================================================================== //
  // Momentum grid
  Grid grid(nfft, idim, dp);
  
  // Green functions
  GF green(grid, ntimes);

  // Self-energies
  Sigma sig(ntimes, grid);

  // Collision terms
  Collision coll(dt, ntimes, grid.nq);

  // Time-stepping
  Propagate td(dt, grid);

  // ======================================================================== //
  //                                  TIME LOOP                               //
  // ======================================================================== //
  int nq = grid.nq;
  double density = 0.0;
  double kinetic = 0.0;
  double amee = 0.067;
  double h2m2 = h2m02 / amee; 
  double scal = 2.0*dp*dp*dp*oneovertpic*aindb*aindb*aindb;
  for(int iq = 0; iq < nq; iq++)
    density += green.gle(iq, 0, 0).imag();
  cout << "density = " << density * scal << endl; 
  for(int iq = 0; iq < nq; iq++)
    kinetic += green.gle(iq, 0, 0).imag() * grid.p2[grid.idxmap[iq]];
  cout << "kinetic = " << kinetic * scal * h2m2 / eryd<< endl; 
  
  sig.calSigfk(grid, green, 0);
  double Pfk = 0.0;
  for(int iq = 0; iq < nq; iq++)
    Pfk += green.gle(iq, 0, 0).imag() * sig.sfk[iq];
  cout << "pe_mf = " << Pfk * scal * 0.5 / eryd<< endl; 

  coll.calColl(green, sig, 0);
  cout << "Igr = " << endl;
  for(int iq = 0; iq < 100; iq++)
    cout << coll.Igr[iq] << "    ";
  cout << endl;
  double Pcor = 0.0;
  for(int iq = 0; iq < nq; iq++)
    Pcor += coll.Ile[iq].imag();
  cout << "pe_cor = " << Pcor * scal * 0.5 / eryd<< endl; 
  
  for(int it = 0; it < ntimes - 1 ; it++){
    
    printf("it = %d\n", it);
    
    // update self-energies
    sig.calSigfk(grid, green, it);
    sig.calSig(grid, green, it);
    cout << "Sigma 1 = " << endl;
    for(int iq = 0; iq < 100; iq++){
      cout << sig.sle[it*nq + iq] << "    ";
    }
    cout << endl;
    cout << "Sigma 2 = " << endl;
    for(int iq = 0; iq < 100; iq++){
      cout << sig.sgr[it*nq + iq] << "    ";
    }
    cout << endl;
    
    Pfk = 0.0;
    for(int iq = 0; iq < nq; iq++)
      Pfk += green.gle(iq, it, it).imag() * sig.sfk[iq];
    cout << "pe_mf = " << Pfk * scal * 0.5 / eryd<< endl; 

    // update collision terms
    coll.calColl(green, sig, it); 
    cout << "Igr = " << endl;
    for(int iq = 0; iq < 100; iq++)
      cout << (coll.Igr[it*nq + iq] - coll.Ile[it*nq + iq]).real() << "    ";
    cout << endl;
    
    Pcor = 0.0;
    for(int iq = 0; iq < nq; iq++)
      Pcor += coll.Ile[it*nq + iq].imag();
    cout << "pe_cor = " << Pcor * scal * 0.5 / eryd<< endl; 
   
    // update Green functions
    td.stepping(green, coll, sig, grid, it);
    
    // output
    density = 0.0;
    for(int iq = 0; iq < nq; iq++)
      density += green.gle(iq, it+1, it+1).imag();
    cout << "density = " << density * scal << endl; 
    
    kinetic = 0.0;
    for(int iq = 0; iq < nq; iq++)
      kinetic += green.gle(iq, it+1, it+1).imag() * grid.p2[grid.idxmap[iq]];
    cout << "kinetic = " << kinetic * scal * h2m2 / eryd<< endl; 
  }

  int n = grid.n;
  complex<double> *fftbox = new complex<double>[n]();
  green.put_into_fftbox(fftbox, grid.idxmap, ntimes-1, ntimes-1);
  for(int i = 0; i < grid.nfft; i++)
    cout << fftbox[i*nfft*nfft].imag() << endl;
  delete[] fftbox;
  
  return 0;
}
