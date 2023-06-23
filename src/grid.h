#ifndef grid_h
#define grid_h

class Grid
{
  public:
    int nfft; // size of the FFT box in each dimension
    int n; // nfft^3
    int idim; // i*i + j*j + k*k <= idim*idim
    int nq; // number of momentum grid
    int *idxmap=nullptr; // nq 
    double dp; // momentum grid interval
               // in Ang^{-1}
    double *mesh=nullptr; // coordinates of momentum grid; nq x 3
                          // in Ang^{-1}
    double *p2 = nullptr; // n
    
    Grid(int nfft_, int idim_, double dp_);
    ~Grid();
};

#endif
