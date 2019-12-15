#ifndef _QUADRATURE_LOBATTO_H_
#define _QUADRATURE_LOBATTO_H_

#include "Quadrature.H"

class Quadrature_Lobatto:public Quadrature
{
  public:
    ///Lobatto for real part of the grid.
    Quadrature_Lobatto(int& a_dvr);
    virtual ~Quadrature_Lobatto();
    virtual std::vector<double> getPoints() const;
    virtual std::vector<double> getWeights() const;
    virtual std::vector<std::vector<double> > getDerivatives() const;
    virtual void print() const;
  private:
    int m_kind, m_kpts, m_dvr;
    double m_alpha, m_beta;
    std::vector<double> m_points, m_weights;
    std::vector<std::vector<double> > m_derivatives;
};
#endif
