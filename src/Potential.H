#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

class Potential
{
public:
  void plot_potential() const;
  std::vector<double> get_potential() const;
protected:
  std::vector<double> m_potential;
  std::vector<double> m_points;
};

class Square: public Potential
{
public:
  Square() {}
  Square(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha, double& V_o);
};

class Step: public Potential
{
public:
  Step() {}
  Step(std::vector<double> i_points, std::vector<double> m_weights, double& V_o);
};

class Delta: public Potential
{
public:
  Delta() {}
  Delta(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha);
};

class Morse: public Potential
{
public:
  Morse() {}
  Morse(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha, double& a_De, double& a_re);
};

class Yukawa: public Potential
{
public:
  Yukawa() {}
  Yukawa(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha);
};

class Coulomb: public Potential
{
public:
  Coulomb() {}
  Coulomb(std::vector<double> i_points, std::vector<double> m_weights, int& a_Z);
};

#endif