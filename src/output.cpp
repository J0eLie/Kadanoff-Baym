#include "output.h"

double Output::write_header()
{
  time_t now = time(NULL);
  cout << "#Time propagation of the Kadanoff-Baym equations begins at "
       << ctime(&now);
  cout << "#Bohr radius (a_B) of exciton in Ang = " << aindb << endl;
  cout << "#One exciton Rydberg (Ry) in eV = " << eryd << endl;
  cout << endl;
  cout << "#Time    Density    Q.Moment    T.Energy    K.Energy    P.Energy (MF)    P.Energy (Corr.)" << endl;
  cout << "#(fs)    (a_B^-3)   (a_B^-2)    (Ry/a_B^3)  (Ry/a_B^3)  (Ry/a_B^3)       (Ry/a_B^3)" << endl; 
  double start_time = omp_get_wtime();
  return start_time;
}

void Output::write_line(int it, const Grid &grid,  const GF &green, const Sigma &sig, const Collision &coll)
{
  t = it*dt;
  
  int nq = grid.nq;
  double dp = grid.dp;
  double h2m2 = h2m02 / amee; 
  double scal = 2.0*dp*dp*dp*oneovertpic*aindb*aindb*aindb;
  
  // density 2/(2*pi)^3 int d^3p f(p, t)
  // f(p, t) = -i*G^<(p, t, t)
  density = 0.0;
  for(int iq = 0; iq < nq; iq++)
    density += green.gle(iq, it, it).imag();
  density *= scal;

  // kinetic energy density 2/(2*pi)^3 int d^3p f(p, t)*h2m2*p2
  KEn = 0.0;
  for(int iq = 0; iq < nq; iq++)
    KEn += green.gle(iq, it, it).imag() * grid.p2[grid.idxmap[iq]];
  KEn *= scal * h2m2 / eryd;

  // mean field(fock) potential energy density 1/(2*pi)^3 int d^3p f(p, t)* Sigma^MF(p, t)
  PEn_mf = 0.0;
  for(int iq = 0; iq < nq; iq++)
    PEn_mf += green.gle(iq, it, it).imag() * sig.sfk[iq];
  PEn_mf *= scal * 0.5 / eryd;

  // correlation potential energy density -i/(2*pi)^3 int d^3p I_1^<(p, t, t)
  PEn_corr = 0.0;
  for(int iq = 0; iq < nq; iq++)
    PEn_corr += coll.Ile[it*nq + iq].imag();
  PEn_corr *= scal * 0.5 / eryd;

  // total energy density
  TEn = KEn + PEn_mf + PEn_corr;

  // quadrupole moment of the distribution
  // Q(t) = 1/density * sqrt(5/16/pi) * 2/(2*pi)^3 int d^3p f(p, t) * [2*pz^2 - px^2 - py^2]
  Qm = 0.0;
  for(int iq = 0; iq < nq; iq++)
    Qm += green.gle(iq, it, it).imag() * (2.0*grid.mesh[iq*3 + 2]*grid.mesh[iq*3 + 2]
                                   - grid.mesh[iq*3]*grid.mesh[iq*3] - grid.mesh[iq*3 + 1]*grid.mesh[iq*3 + 1]);
  Qm *= scal * aindb * aindb * sqrt(5.0/16.0/pi) / density;

  // write to stdout
  cout << setiosflags(ios::right) << setprecision(6) << setiosflags(ios::fixed);
  cout << setw(12) << t << setw(12) << density << setw(12) << Qm << setw(12)
       << setw(12) << TEn << setw(12) << KEn << setw(12) << PEn_mf << setw(12) << PEn_corr << endl;

}

void Output::write_distrib(int it, const Grid &grid, const GF &green)
{
  if(it % outfreq == 0){
    double *fp = new double[grid.n]();
    for(int iq = 0; iq < grid.nq; iq++)
      fp[grid.idxmap[iq]] = green.gle(iq, it, it).imag();
    ofstream of("distribution_" + to_string(it), ios::out);
    int nfft = grid.nfft;
    of << setiosflags(ios::right) << setprecision(6) << setiosflags(ios::fixed);
    for(int i = 0; i < nfft; i++){
      for(int j = 0; j < nfft; j++){
        for(int k = 0; k < nfft; k++) of << setw(10) << fp[i*nfft*nfft + j*nfft + k];
      }
    }
    of.close();

    delete[] fp;
  }
}

void Output::write_ender(double start_time)
{
  time_t now = time(NULL);
  cout << endl;
  cout << "#End of propagation at " << ctime(&now);
  double end_time = omp_get_wtime();
  cout << "#Total time used " << end_time - start_time << " seconds" << endl;
}

Output::Output(double dt_, int pfreq)
{
  dt = dt_;
  outfreq = pfreq;
}
Output::~Output() {}
