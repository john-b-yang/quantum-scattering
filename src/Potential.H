#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

/// The base class of Potential
class Potential
{
public:
  /// Plot function of Potential Class
  void plot_potential() const;
  /// Get function of Potential Class
  std::vector<double> get_potential() const;
protected:
  std::vector<double> m_potential;
  std::vector<double> m_points;
};

/// A Square class derived from the base Potential class
class Square: public Potential
{
public:
  /// Default constructor
  Square() {}
  /// Constructor of Square
  Square(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha, double& V_o);
};

/// A Step class derived from the base Potential class
class Step: public Potential
{
public:
  /// Default constructor
  Step() {}
  /// Constructor of Step
  Step(std::vector<double> i_points, std::vector<double> m_weights, double& V_o);
};

/// A Delta class derived from the base Potential class
class Delta: public Potential
{
public:
  /// Default constructor
  Delta() {}
  /// Constructor of Delta
  Delta(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha);
};

/// A Morse class derived from the base Potential class
class Morse: public Potential
{
public:
  /// Default constructor
  Morse() {}
  /// Constructor of Morse
  Morse(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha, double& a_De, double& a_re);
};

/// A Yukawa class derived from the base Potential class
class Yukawa: public Potential
{
public:
  /// Default constructor
  Yukawa() {}
  /// Constructor of Yukawa
  Yukawa(std::vector<double> i_points, std::vector<double> m_weights, double& a_alpha);
};

/// A Coulomb class derived from the base Potential class
class Coulomb: public Potential
{
public:
  /// Default constructor
  Coulomb() {}
  /// Constructor of Coulomb
  Coulomb(std::vector<double> i_points, std::vector<double> m_weights, int& a_Z);
};

#endif