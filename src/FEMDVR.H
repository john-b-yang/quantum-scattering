#ifndef _FEMDVR_H_
#define _FEMDVR_H_
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>

#include "Quadrature_Lobatto.H"
#include "Quadrature_Radau.H"

class FEMDVR
{
  public:
    /// Interface base class for FEM-DVR grid 

    /// Default constructor 
    FEMDVR();

    /// Constructor arguments: lobatto and radau order, number of elements.
    FEMDVR(int& a_Lndvr, int& a_Rndvr, int& a_Nelem)
      : m_theta(30.0), m_alphaRad(0.3)
    {
      m_Lndvr = a_Lndvr;
      m_Rndvr = a_Rndvr;
      m_Nelem = a_Nelem;
      m_Nbas = (a_Nelem-1)*(a_Lndvr-1)+a_Rndvr-1;
      m_eit = std::polar(1.0,M_PI*m_theta/180.0);

      //Default finite element spacing
      for(int i = 0; i < a_Nelem+1; ++i)
        {
          m_realbounds.push_back(i*10);
        }
      m_R0 = m_realbounds[-1];
      m_complexbounds = m_realbounds;
      m_complexbounds.back() = m_R0 + m_eit*(m_realbounds.back()-m_R0);

      Quadrature_Lobatto lobattoQuad(a_Lndvr);
      std::vector<double> lobX = lobattoQuad.getPoints();
      std::vector<double> lobW = lobattoQuad.getWeights();
      for (int element = 0; element < m_Nelem-1; ++element)
        {
          std::vector<std::complex<double> > lobtempX, lobtempW;
          for (int dvr = 0; dvr < lobX.size(); ++dvr)
            {
              lobtempX.push_back(m_complexbounds[element]+0.5*(lobX[dvr]+1)*
                      (m_complexbounds[element+1]-m_complexbounds[element]));
              lobtempW.push_back(0.5*lobW[dvr]*
                      (m_complexbounds[element+1]-m_complexbounds[element]));
            }
          if (element!=0) m_weights.back()=m_weights.back()+lobtempW.front();
          for (int dvr = 1; dvr < lobtempX.size(); ++dvr)
            {
              m_points.push_back(lobtempX[dvr]);
              m_weights.push_back(lobtempW[dvr]);
            }
        }

      Quadrature_Radau radauQuad(a_Rndvr); 
      std::vector<double> radX = radauQuad.getPoints();
      std::vector<double> radW = radauQuad.getWeights();
      std::vector<std::complex<double> > radtempX, radtempW;
      std::complex<double> radscale = m_eit/(2.0*m_alphaRad);
      for (int i = 0; i < radX.size(); ++i)
        {
          radtempX.push_back(radscale*radX[i] + m_R0);
          radtempW.push_back(radscale*radW[i]);
        }
      m_weights.back()=m_weights.back()+radtempW.front();
      for (int i = 1; i < radtempX.size(); ++i)
        {
          m_points.push_back(radtempX[i]);
          m_weights.push_back(radtempW[i]);
        }

      m_complexbounds[m_Nelem] = radtempX[a_Rndvr];
      m_realbounds[m_Nelem] = (m_complexbounds[m_Nelem]-m_R0)*
                              conj(m_eit)+m_R0;
    } //end of constructor

    /// Destructor.
    /* virtual ~FEMDVR(){}; */

    /// Getter for points.
    std::vector<std::complex<double> > getPoints() const
    {
      return m_points;
    }

    /// Getter for weights.
    std::vector<std::complex<double> > getWeights() const
    {
      return m_weights;
    }

    /// Plots the DVR basis functions in DIM dimensions.

    /* void plotDVR(int DIM) const; */

    /// Get number of basis functions
    const int& getNbas(){return m_Nbas;}

    /// Get R0
    const std::complex<double>& getR0(){return m_R0;}

    virtual void print() const
    {
      std::cout << std::endl;
      std::cout << "*** FEM-DVR Grid***" << std::endl;
      std::cout << std::endl;
      std::cout << "Lobatto Order = " << m_Lndvr << "; " 
                << "Radau Order = " << m_Rndvr << std::endl;
      std::cout << "Number of Finite-Elements = " << m_Nelem << "; "
                << "Number of Basis Functions = " << m_Nbas << std::endl;
      std::cout << std::endl;
      for( int i = 0; i<m_points.size(); ++i )
        {
          std::cout << "Point " << m_points[i] << "; " 
                    << "Weights " << m_weights[i] << std::endl;
        }
      std::cout << "R0 = " << m_R0 << "; "
                << "Theta = " << m_theta << std::endl;
    }

  protected:
    int m_Nelem, m_Nbas, m_Lndvr, m_Rndvr;
    double m_theta, m_alphaRad;
    std::complex<double> m_eit, m_R0;
    std::vector<std::complex<double> > m_points, m_weights;
    std::vector<std::complex<double> > m_realbounds, m_complexbounds;
};
#endif
