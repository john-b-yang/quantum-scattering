#ifndef _QUADRATURE_LOBATTO_H_
#define _QUADRATURE_LOBATTO_H_

extern "C" 
{
  void gaussq_(int*, int*, double*, 
              double*, int*, double*, 
              double*, double*, double*);
}

class Quadrature_Lobatto
{
  public:
    ///Lobatto wrapper to Fortran quadrature routines
    Quadrature_Lobatto();
    Quadrature_Lobatto(int a_dvr);
    ~Quadrature_Lobatto();
    std::vector<double> getPoints();
    std::vector<double> getWeights();
    void print() const;
  private:
    unsigned int m_dvr;
    int m_kind, m_kpts;
    double m_alpha, m_beta;
    std::vector<double> m_points, m_weights;
};
#endif
