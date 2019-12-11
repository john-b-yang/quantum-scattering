#include <vector>
#include <cmath>
#include <iostream>
#include <stdio.h>

#include "Quadrature_Lobatto.H"

Quadrature_Lobatto::Quadrature_Lobatto(){}

Quadrature_Lobatto::Quadrature_Lobatto(int a_dvr)
  : m_kind(1), m_kpts(2), m_alpha(0), m_beta(0)

{
  m_dvr = a_dvr;
  double* x = new double[a_dvr];
  double* w = new double[a_dvr];
  double* src = new double[a_dvr];
  double* endpts = new double[2];
  endpts[0] = -1.0; endpts[1] = 1.0;
  gaussq_(&m_kind, &a_dvr, &m_alpha, &m_beta, &m_kpts, endpts, 
          src, x, w);
  for (int i=0; i<a_dvr; ++i)
    {
      std::cout << x[i] << std::endl;
      m_points.push_back(x[i]);
      m_weights.push_back(w[i]);
    }

}

std::vector<double> Quadrature_Lobatto::getPoints()
{
  return m_points;
}

std::vector<double> Quadrature_Lobatto::getWeights()
{
  return m_weights;
}

void Quadrature_Lobatto::print() const
{
  for(int i=0; i<m_dvr; ++i)
    {
      std::cout<< "Points " << m_points[i] << " " << " Weights" << m_weights[i] << std::endl;
    }
}
