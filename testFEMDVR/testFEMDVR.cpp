#include <iostream>
#include <vector>

#include "FEMDVR.H"

int main(int argc, char** argv)
{
  int ndvr = 5, mdvr = 5, nel = 3;

  FEMDVR FEMDVR(ndvr, mdvr, nel);

  FEMDVR.print();
  /* FEMDVR.plotDVR(); */

  return 0;
}


