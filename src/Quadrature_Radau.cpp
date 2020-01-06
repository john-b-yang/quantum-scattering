#include <vector>
#include <cmath>
#include <iostream>
#include <stdio.h>

#include "Quadrature.H"
#include "Quadrature_Radau.H"

Quadrature_Radau::Quadrature_Radau(int& a_dvr):Quadrature(a_dvr)
{
  m_kind = 6, m_kpts = 1, m_alpha = 0, m_beta = 0;
  double* x = new double[a_dvr];
  double* w = new double[a_dvr];
  double* src = new double[a_dvr];
  double* endpts = new double[2];
  endpts[0] = 0.0; endpts[1] = 1.0;
  gaussq_(&m_kind, &a_dvr, &m_alpha, &m_beta, &m_kpts, endpts, 
          src, x, w);
  for (int i=0; i<a_dvr; ++i)
    {
      m_points.push_back(x[i]);
      m_weights.push_back(w[i]);
    }
  m_derivatives.resize(a_dvr);
  for (int i=0; i<a_dvr; ++i)
    {
      m_derivatives[i].resize(a_dvr);
      m_derivatives[i][i] = 0.0;
      for (int k=0; k<a_dvr; ++k)
        {
          if(i!=k) m_derivatives[i][i] += 1.0/(x[i]-x[k]);
        }
      for (int j=0; j<a_dvr; ++j)
        {
          if(i!=j)
            {
              double temp = 1.0/(x[j]-x[i]);
              for (int k=0; k<a_dvr; ++k)
                {
                  if((k!=i)&&(k!=j)) temp *= (x[i]-x[k])/(x[j]-x[k]);
                }
              m_derivatives[i][j] = temp;
            }
        }
    }
}

Quadrature_Radau::~Quadrature_Radau(){};

std::vector<double> Quadrature_Radau::getPoints() const
{
  return m_points;
}

std::vector<double> Quadrature_Radau::getWeights() const
{
  return m_weights;
}

std::vector<std::vector<double> > Quadrature_Radau::getDerivatives() const
{
  return m_derivatives;
}

void Quadrature_Radau::print() const
{
  std::cout<< "Radau Quadrature" << std::endl;
  for(int i=0; i<m_dvr; ++i)
    {
      std::cout<< "Points " << m_points[i] << " " << 
                  "Weights " << m_weights[i] << std::endl;
    }
  for(int i=0; i<m_dvr; ++i)
    {
      for(int j=0; j<m_dvr; ++j)
        {
          std::cout << "Derivates " << m_derivatives[i][j] << std::endl;
        }
    }
}
