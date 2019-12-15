#ifndef _QUADRATURE_RADAU_H_
#define _QUADRATURE_RADAU_H_

#include "Quadrature.H"

class Quadrature_Radau:public Quadrature
{
  public:
    ///Radau for complex tail. 
    Quadrature_Radau(int& a_dvr);
    virtual ~Quadrature_Radau();
    virtual std::vector<double> getPoints() const;
    virtual std::vector<double> getWeights() const;
    virtual std::vector<std::vector<double> > getDerivatives() const;
    virtual void print() const;
  private:
    int m_kind, m_kpts;
    double m_alpha, m_beta;
    std::vector<double> m_points, m_weights;
    std::vector<std::vector<double> > m_derivatives;
};
#endif
