#ifndef _FEMDVR_H_
#define _FEMDVR_H_

#include "Potential.H"
class StaticPotential

class FEMDVR
{
  public:
    /// Sets FEM-DVR grid with Hamiltonian representation for a given potential.  Also has PR_TIMERS integrated.

    /// Default constructor 
    FEMDVR();

    /// Class that constructs FEM-DVR grid and grid representation of Hamiltonian
    FEMDVR(int a_nel, const StaticPotential& a_potential);
      : m_theta(30.0), m_alphaRad(0.3)

    /// Destructor
    ~FEMDVR();

    /// Getter for points
    vector<double> getPoints() const;

    /// Getter for weights
    vector<double> getWeights() const;

    /// Plots the DVR basis functions in DIM dimensions.
    void plotDVR(int DIM) const;

    /// Plots the Wavefunction on the FEM-DVR grid in DIM dimensions.
    void plotWave(int DIM) const;

    /// Checks to see whether the Hamiltonian is symmetric.
    bool symmetric() const;

    /// Get number of basis functions
    int getNbas() const;

    /// Get Hamiltonian
    double getTXX() const;

    /// Get Inverse Hamiltonian
    double getTIXX() const;

    /// Get R0
    double getR0() const;

    void print() const;
  private:
    unsigned int m_nel, m_nbas, lmax;
    double m_theta, m_alphaRad, m_R0;
    complex m_eit;
    vector<double> m_realBoundaries;
    vector<complex> m_complexBoundaries;
    vector<vector<vector<double> > > TXX, TIXX; // Hamiltonian and inverse
}
#endif
