#include <cstdio>   
#include <cmath>    
#include <cstring>  
#include <cassert>  
#include <vector>

#include "Quadrature_Lobatto.H"
/* #include "Quadrature_Radau.H" */

int main(int argc, char** argv)
{
  int ndvr = 5, mdvr = 5;
  
  Quadrature_Lobatto gridLobatto(ndvr);
  /* Quadrature_Radau gridRadau(mdvr); */
 
  /* gridLobatto.print(); */
  /* gridRadau.print(); */

  return 0;
}

  
