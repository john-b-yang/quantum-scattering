#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

extern "C" 
{
  void gaussq_(int*, int*, double*, 
              double*, int*, double*, 
              double*, double*, double*);
}

class Quadrature
{
  public:
    ///Quadrature base class that wraps fortran routine gaussq
    Quadrature(){};
    Quadrature(int a_dvr){m_dvr = a_dvr;}
    virtual ~Quadrature(){};
    virtual std::vector<double> getPoints() const = 0;
    virtual std::vector<double> getWeights() const = 0;
    virtual std::vector<std::vector<double> > getDerivatives() const = 0;
    virtual void print() const = 0;

    //Access function
    const unsigned int& getDVROrder(){return m_dvr;}
  protected:
    unsigned int m_dvr;
};
#endif
