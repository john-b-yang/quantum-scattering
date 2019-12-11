#include <iostream>   
#include <vector>

#include "FEMDVR.H"
#include "StaticPotential.H"

int main(int argc, char** argv)
{
  int ndvr = 5, mdvr = 5, nel = 2;
  double r;

  r = 5.0;
  StaticPotential.square(r);
  
  FEMDVR FEMDVR(nel, potential);
 
  FEMDVR.print();

  return 0;
}

  
