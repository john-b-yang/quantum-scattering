#ifndef _COUL90_H_
#define _COUL90_H_

#include <vector>
#include <complex>

// We maynot use this but it would be cool. Can at least put it on the DVR grid.

class COUL90
{
  public:
    /// A wrapper for the COUL90 Fortran routines that computes hypergeometric functions
    COUL90();

    COUL90(int a_nel);// need to figure out the correct intializing parameters

    /// Generates the expansion of a coulomb wave on the DVR grid
    void coulombDVR() const;

    /// Generates the one-electron testing functions expansion coefficients on DVR grid
    void orbDVR() const;

    void plotCoulomb() const;

  private:
    int m_info, m_lmax;
    double m_k, m_reta;
    vector<double> m_fc, m_fcp, m_gc, m_gcp
    vector<vector<complex<double> > > m_coulexp, m_psik;
}
#endif
