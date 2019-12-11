# Fortran code
  * gq.f90 is what needs to be imported to get Lobatto quadratures
  * quads.f90 shows an interface to gq.f90 to get different quadratures
  * build_grids.f90 builds the FEM-DVR grid
  * cdiag.f90 diagonalizes a matrix to give eigenvalues and eigenvectors
  * COUL90.f is a fortran 77 code that computes a coulomb wavefunction
  * build_KEmats.f90 builds the Hamiltonian on the DVR grid
  * build_R12matsRad.f90 builds the two-electron $\frac{1}{\abs{r_1-r_2}}$ term

# Laugerre PDF
  * This PDF explains the Radau tail in detail
  * All we care about is how to transform our Radau raw points to the FEM grid
    which is found in equation $32$ on page $6$.
  
