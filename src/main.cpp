#include <string>
using namespace std;

#include "grid.h"
#include "gf.h"
#include "sigma.h"
#include "collision.h"
#include "propagate.h"
#include "output.h"
#include "constants.h"

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
  double dp = 0.3/aindb;
  
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

  // Output
  Output out(dt, 15);

  // ======================================================================== //
  //                                  TIME LOOP                               //
  // ======================================================================== // 
  
  double start_time = out.write_header();

  for(int it = 0; it < ntimes - 1 ; it++){
    
    // update self-energies
    sig.calSigfk(grid, green, it);
    sig.calSig(grid, green, it);
    
    // update collision terms
    coll.calColl(green, sig, it); 
    
    if(it == 0) out.write_line(0, grid, green, sig, coll);

    // update Green functions
    td.stepping(green, coll, sig, grid, it);
    
    // output
    out.write_line(it + 1, grid, green, sig, coll);
  }
  
  out.write_ender(start_time);

  return 0;
}
