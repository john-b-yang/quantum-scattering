#ifndef _QUADRATURE_RADAU_H_
#define _QUADRATURE_RADAU_H_

extern "C"  
{
  void gaussq_(int*, int*, double*, 
              double*, int*, double*, 
              double*, double*, double*);
}

class Quadrature_Radau
{
  public:
    ///Radau wrapper to Fortran quadrature routines
    Quadrature_Radau();
    Quadrature_Radau(int a_dvr);
    ~Quadrature_Radau();
    std::vector<double> getPoints() const;
    std::vector<double> getWeights() const;
    void print() const;
  private:
    unsigned int m_dvr;
    int m_kind, m_kpts;
    double m_alpha, m_beta;
    std::vector<double> m_points, m_weights;
};
#endif
