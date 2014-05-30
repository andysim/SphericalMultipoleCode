Tinker format)# Spherical Multipole Code

This repository provides a set of Fortran routines that aid in the implementation of the algorithms detailed in [our recent JCP article](http://scitation.aip.org/content/aip/journal/jcp/140/18/10.1063/1.4873920).  The analytic form of all equations can be found in the [supplementary material](ftp://ftp.aip.org/epaps/journ_chem_phys/E-JCPSA6-140-005418) for that paper.  To aid in the debugging of these routines' implementation, we provide detailed output showing the intermediates involved for a hexadecapole-bearing water dimer, with the following parameters.

##Geometry (in Angstroms)

         O    2.000000    2.000000    2.000000
         H    2.500000    2.000000    3.000000
         H    1.500000    2.000000    3.000000
         O    0.000000    0.000000    0.000000
         H    0.500000    0.000000    1.000000
         H   -0.500000    0.000000    1.000000

##Multipoles on each center

The ordering and normalization follows the Tinker convention, which uses atomic units.  The explicit order is
Dipoles       x    y    z
Quadrupoles   xx   xy   yy   xz   yz   zz
Octopoles     xxx  xxy  xyy  yyy  xxz  xyz  yyz  xzz  yzz  zzz
Hexadecapoles xxxx xxxy xxyy xyyy yyyy xxxz xxyz xyyz yyyz xxzz xyzz yyzz xzzz yzzz zzzz

Our toy system places no multipoles on the hydrogens, but uses the following hexadecapoles on oxygen, oriented using the bisector of the hydrogen atoms.

Dipoles          0.0000  0.0000  0.7283
Quadrupoles      0.2432  0.0000 -0.1350  0.0000  0.0000  -0.1082
Octopoles        0.0000  0.0000  0.0000  0.0000  3.2459   0.0000  -1.3462  0.0000  0.0000  -1.8997
Hexadecapoles   -3.4510  0.0000 -0.7318  0.0000  1.4036   0.0000   0.0000  0.0000  0.0000   4.1829  0.0000  -0.6718  0.0000  0.0000  -3.5111


