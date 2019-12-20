#include <fstream>
#include <vector>
#include <iostream>
#include <Potential.H>
#include <cmath>
#define Const_E        2.71828182845904523536

void Potential::plot_potential() const
{
  system("mkdir -p Data");
  std::ofstream outfile;
  outfile.open("potential.dat");
  for (int i = 0; i < m_potential.size(); ++i)
  {
    outfile << m_points[i] << " " 
                  << m_potential[i] << std::endl;
  }
  outfile.close();
  system("mv potential.dat Data/");

  system("mkdir -p Figures");
  gnuplot plot;
  plot("set terminal postscript eps color");
  plot("set output \"Potential.eps\" ");
  plot("plot \"./Data/potential.dat\" u 1:2 w l ");
  system("mv Potential.eps Figures/");
}

std::vector<double> Potential::get_potential() const
{
  return m_potential;
}

Square::Square(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha, double& V_o)
{
  for (int i = 0; i < i_points.size(); i++) {
    double tem = i_points[i], resp = 0;
    m_points.push_back(tem);
    if (tem < 0 || tem > a_alpha) {
      resp = 0;
    }
    else {
      resp = V_o;
    }
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Step::Step(std::vector<double> i_points, std::vector<double> m_weights, double& V_o)
{
  for (int i = 0; i < i_points.size(); i++) {
    double tem = i_points[i], resp = 0;
    m_points.push_back(tem);
    if (tem < 0) {
      resp = 0;
    }
    else {
      resp = V_o;
    }
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Delta::Delta(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha)
{
  for (int i = 0; i < i_points.size(); i++) {
    double tem = i_points[i], resp = 0;
    m_points.push_back(tem);
    if (tem == a_alpha) {
      resp = 1;
    }
    else {
      resp = 0;
    }
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Morse::Morse(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha, double& a_De, double& a_re)
{
  for (int i = 0; i < i_points.size(); i++) {
    double tem = i_points[i], resp = 0;
    m_points.push_back(tem);
    double tem2 = 1 - std::pow(Const_E, -a_alpha*(tem - tem * a_re));
    resp = a_De * tem2 * tem2;
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Yukawa::Yukawa(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha)
{
  for (int i = 0; i < i_points.size(); i++) {
    double tem = i_points[i], resp = 0;
    m_points.push_back(tem);
    resp = std::pow(Const_E, -(a_alpha * tem)) / tem;
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Coulomb::Coulomb(std::vector<double> i_points, std::vector<double> m_weights, int& a_Z)
{
  for (int i = 0; i < i_points.size(); i++) {
    double tem = i_points[i], resp = 0;
    m_points.push_back(tem);
    resp = a_Z / tem;
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}