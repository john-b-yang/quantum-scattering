#include <vector>
#include <iostream>
#include <Potential.H>
#include <cmath>
#define Const_E        2.71828182845904523536

void Potential::plot_potential() const
{

}

std::vector<double> Potential::get_potential() const
{
  return m_potential;
}

Square::Square(std::vector<double> m_points, std::vector<double> m_weights, double& a_r)
{
  for (int i = 0; i < m_points.size(); i++) {
    double tem = m_points[i], resp = 0;
    if (tem < 0) {
      resp = 0;
    }
    else {
      //resp = V_o; ??
    }
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Step::Step(std::vector<double> m_points, std::vector<double> m_weights, double& a_r)
{
  for (int i = 0; i < m_points.size(); i++) {
    double tem = m_points[i], resp = 0;
    if (tem < 0) {
      resp = 0;
    }
    else {
      //resp = V_o; ??
    }
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Delta::Delta(std::vector<double> m_points, std::vector<double> m_weights, double& a_alpha, double& a_r)
{
  for (int i = 0; i < m_points.size(); i++) {
    double tem = m_points[i], resp = 0;
    resp = a_alpha * (a_r);
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Morse::Morse(std::vector<double> m_points, std::vector<double> m_weights, double& a_alpha, double& a_De, double& a_re, double& a_r)
{
  for (int i = 0; i < m_points.size(); i++) {
    double tem = m_points[i], resp = 0;
    double tem2 = 1 - std::pow(Const_E, -a_alpha*(a_r - a_re));
    resp = a_De * tem2 * tem2;
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Yukawa::Yukawa(std::vector<double> m_points, std::vector<double> m_weights, double& a_alpha, double& a_r)
{
  for (int i = 0; i < m_points.size(); i++) {
    double tem = m_points[i], resp = 0;
    resp = std::pow(Const_E, -(a_alpha * a_r)) / a_r;
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}

Coulomb::Coulomb(std::vector<double> m_points, std::vector<double> m_weights, int& a_Z, double& a_r)
{
  for (int i = 0; i < m_points.size(); i++) {
    double tem = m_points[i], resp = 0;
    resp = a_Z / a_r;
    resp *= m_weights[i];
    m_potential.push_back(resp);
  }
}