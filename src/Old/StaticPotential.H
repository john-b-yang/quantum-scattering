#ifndef _STATICPOTENTIAL_H_
#define _STATICPOTENTIAL_H_

class StaticPotential
{
  public:

    StaticPotential();

    ///Square well static Potential $=0$ if $x<0, \alpha<x$ and $=V_o$ if $0<x\alpha$
    void square(double& a_r);

    ///Step static Potential $=0$ if $x<0$ and $=V_o$ if $x\geq0$
    void step(double& a_r);

    ///Delta static Potential $\alpha \delta(r)$
    void delta(double& a_alpha, double& a_r);

    ///Morse static Potential $De(1 - e^{-\alpha(r-re)})^2$
    void morse(double& a_alpha, double& a_De, double& a_re, double& a_r);

    ///Yukawa static Potential $\frac{e^{-\alpha r}}{r}$
    void yukawa(double& a_alpha, double& a_r);

    ///Coulomb Potential $\frac{Z}{r}$
    void coulomb(int& a_Z, double& a_r);

  private:
}
#endif
