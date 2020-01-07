#ifndef _FEMDVR_H_
#define _FEMDVR_H_
#include <cmath>
#include <vector>
#include <complex>
#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "Quadrature_Lobatto.H"
#include "Quadrature_Radau.H"
#include "gnuplot.H"

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
      m_NRealbas = (a_Nelem-1)*(a_Lndvr-1)-1;
      m_Nbas = (a_Nelem-1)*(a_Lndvr-1)+a_Rndvr-1;
      m_eit = std::polar(1.0,M_PI*m_theta/180.0);

      //Default finite element spacing
      for(int i = 0; i < a_Nelem+1; ++i)
        {
          m_realbounds.push_back(i*10);
        }
      m_R0 = m_realbounds[a_Nelem-1];
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
      for (int dvr = 1; dvr < radtempX.size(); ++dvr)
        {
          m_points.push_back(radtempX[dvr]);
          m_weights.push_back(radtempW[dvr]);
        }

      m_complexbounds[m_Nelem] = radtempX[a_Rndvr];
      m_realbounds[m_Nelem] = (m_complexbounds[m_Nelem]-m_R0)*
                              conj(m_eit)+m_R0;

    } //end of constructor for FEMDVR grid

    /// Destructor.
    virtual ~FEMDVR(){};

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

    /// Evaluates function on FEM-DVR grid.
    std::complex<double> dvrEvaluate(const double& a_x, const std::vector<std::complex<double> >& a_func)
    {
      assert(a_x>m_complexbounds[0].real() && a_x<=m_complexbounds[m_Nelem].real());
      int FEM_loc = 0;
      for(int element = 0; element < m_Nelem; ++element)
        {
          if(a_x<element)
            {
              ++FEM_loc;
            }
          else break;
        }
      int index = 0, start = 0 , end = 0;
      std::complex<double>fac = 0;
      std::vector<std::complex<double> > tempPoints, tempWeights;
      if(FEM_loc==1)
        {
          index=1; start=1; end=m_Lndvr; fac=1.0;
          tempPoints.push_back(0.0);
          tempWeights.push_back(0.0);
          std::copy(m_points.begin() + 1, m_points.begin() + m_Lndvr - 1,
                    std::back_inserter(tempPoints));
          std::copy(m_weights.begin() + 1, m_weights.begin() + m_Lndvr - 1,
                    std::back_inserter(tempWeights));
        }
      else if ((FEM_loc!=1) && (FEM_loc!=m_Nelem))
        {
          index=(FEM_loc-1)*(m_Lndvr-1); start=0; end=m_Lndvr; fac=1.0;
          std::copy(m_points.begin() + (FEM_loc)*(m_Lndvr-2),
                    m_points.begin() +  (FEM_loc+1)*(m_Lndvr-1),
                    std::back_inserter(tempPoints));
          std::copy(m_weights.begin() + (FEM_loc)*(m_Lndvr-2),
                    m_weights.begin() +  (FEM_loc+1)*(m_Lndvr-1),
                    std::back_inserter(tempWeights));
        }
      else
        {
          index=m_Nbas-m_Rndvr+1; start=0; end=m_Rndvr;
          fac = std::exp(-m_alphaRad*conj(m_eit)*(a_x-m_R0));
          std::copy(m_points.begin() + m_Nbas - m_Rndvr,
                    m_points.begin() + m_Nbas,
                    std::back_inserter(tempPoints));
          std::copy(m_weights.begin() + m_Nbas - m_Rndvr,
                    m_weights.begin() + m_Nbas,
                    std::back_inserter(tempWeights));
        }
      std::complex<double> accum = 0.0;
      for (int j = start; j < end; ++j)
        {
          std::complex<double> prod = 1.0;
          for (int k = start; k < end; ++k)
            {
              if (k != j)
                {
                  prod *= (a_x-tempPoints[k])/
                          (tempPoints[j]-tempPoints[k]);
                }
            }
          accum += prod*a_func[index]/sqrt(tempWeights[index]);
        }
      return accum*fac;
    }

    /// Plots a function in the DVR representation.
    void plotFunction(const std::vector<std::complex<double> >& a_func, const double& a_rmax, const int& a_npoints ) const
    {
      double dr = a_rmax/a_npoints;
      std::complex<double> zrad, fval;
      system("mkdir -p Data");
      std::ofstream outfile;
      outfile.open("dvrFunc.dat");
      for (int ir = 0; ir < a_npoints; ++ir)
        {
          zrad = irs*dr;
          fval = dvrEvaluate(zrad, a_func);
          for (int i = 0; i< m_points.size(); ++i)
            {
              outfile << zrad << " "
                      << fval << std::endl;
            }
        }
      outfile.close();
      system("mv dvrFunc.dat Data/");

      system("mkdir -p Figures");
      gnuplot plot;
      plot("set terminal postscript eps color");
      plot("set output \"DVRFunc.eps\" ");
      plot("plot \"./Data/dvrFunc.dat\" u 1:2 w l ");
      system("mv DVRFunc.eps Figures/");
    }

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
                    << "Weights " << m_weights[i] << "; "
                    << "DVR Basis " << m_dvrBasis[i] << std::endl;
        }
      std::cout << "R0 = " << m_R0 << "; "
                << "Theta = " << m_theta << std::endl;
    }

  protected:
    int m_Nelem, m_Nbas, m_NRealbas, m_Lndvr, m_Rndvr;
    double m_theta, m_alphaRad;
    std::complex<double> m_eit, m_R0;
    std::vector<std::complex<double> > m_points, m_weights;
    std::vector<std::complex<double> > m_realbounds, m_complexbounds;
    std::vector<std::complex<double> > m_dvrBasis;
};
#endif
