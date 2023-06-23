#include "coulomb.h"

void Coulomb::potcoul(Grid &grid) { // V(p) = 4*pi*e2/eps_b / (p^2 + kappa^2) 
  int n = grid.n;
  V = new double[n];
  double coul = fpi*e2/epsilon_b;
  for(int iq = 0; iq < n; iq++){
    V[iq] = coul / (grid.p2[iq] + kappa*kappa);
  }

  return;
}

void Coulomb::potcoul_fk(Grid &grid) { // Vfk(p) = 4*pi*e2/eps_b / p^2 , |p| > 0
                                       //        = 8*pi*e2/eps_b * (pi + 0.69546)/dp^2, |p| = 0
  int n = grid.n;
  Vfk = new double[n];
  double coul = fpi*e2/epsilon_b;
  Vfk[0] = 2.0 * coul * (pi + 0.69546) / (grid.dp*grid.dp);
  for(int iq = 1; iq < n; iq++){
    Vfk[iq] = coul / grid.p2[iq];
  }

  return;
}

void Coulomb::calakp(GF &gf, Grid &grid) {
  return;
}

Coulomb::Coulomb(Grid &grid) {
  potcoul(grid);
  potcoul_fk(grid);
}

Coulomb::~Coulomb() {
  delete[] V;
  delete[] Vfk;
}
