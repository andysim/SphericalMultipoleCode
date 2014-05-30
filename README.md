# Spherical Multipole Code

This repository provides a set of Fortran routines that aid in the implementation of the algorithms detailed in [our recent JCP article](http://scitation.aip.org/content/aip/journal/jcp/140/18/10.1063/1.4873920).  The analytic form of all equations can be found in the [supplementary material](ftp://ftp.aip.org/epaps/journ_chem_phys/E-JCPSA6-140-005418) for that paper.  We will use notation [proposed by Hättig](http://dx.doi.org/10.1016/S0009-2614(97)00206-6) in his original derivation of the unscreened method.  To aid in the debugging of these routines' implementation, we provide detailed output showing the intermediates involved for a hexadecapole-bearing water dimer, with the following parameters.

Geometry (in Angstroms)
=======================

         O    2.000000    2.000000    2.000000
         H    2.500000    2.000000    3.000000
         H    1.500000    2.000000    3.000000
         O    0.000000    0.000000    0.000000
         H    0.500000    0.000000    1.000000
         H   -0.500000    0.000000    1.000000

Multipoles on each center
=========================

The ordering and normalization follows the Tinker convention, which uses atomic units.  The explicit order is

              |   1  |   2  |   3  |   4  |   5  |   6  |   7  |   8  |   9  |  10  |  11  |  12  |  13  |  14  |  15 
--------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|----:
Dipoles       |    x |    y |    z |      |      |      |      |      |      |      |      |      |      |      |     
Quadrupoles   |   xx |   xy |   yy |   xz |   yz |   zz |      |      |      |      |      |      |      |      |     
Octopoles     |  xxx |  xxy |  xyy |  yyy |  xxz |  xyz |  yyz |  xzz |  yzz |  zzz |      |      |      |      |     
Hexadecapoles | xxxx | xxxy | xxyy | xyyy | yyyy | xxxz | xxyz | xyyz | yyyz | xxzz | xyzz | yyzz | xzzz | yzzz | zzzz

Our toy system places no multipoles on the hydrogens, but uses the following hexadecapoles on oxygen, oriented using the bisector of the hydrogen atoms.

              |   1  |   2  |   3  |   4  |   5  |   6  |   7  |   8  |   9  |  10  |  11  |  12  |  13  |  14  |  15 
--------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|----:
Dipoles       |  0.0000 | 0.0000 |  0.7283
Quadrupoles   |  0.2432 | 0.0000 | -0.1350 | 0.0000 | 0.0000 | -0.1082
Octopoles     |  0.0000 | 0.0000 |  0.0000 | 0.0000 | 3.2459 |  0.0000 | -1.3462 | 0.0000 | 0.0000 | -1.8997
Hexadecapoles | -3.4510 | 0.0000 | -0.7318 | 0.0000 | 1.4036 |  0.0000 |  0.0000 | 0.0000 | 0.0000 |  4.1829 | 0.0000 | -0.6718 | 0.0000 | 0.0000 | -3.5111

The attenuation parameter, kappa, is 0.37 A^{-1}


## Procedure

 - Parse Cartesian moments, convert to appropriate units (we use AKMA), and take linear combinations to form sphericals.
 - Before the pair loops, find the local orientation matrix _Ua_ that orients each center, _a_.
 - For each center, use the _U_ matrix (which rotates dipoles) to build rotation matrices for the higher-order multipoles.
 - Use these rotation matrices to rotate each multipole to the lab frame, and store for later use.
 - For each pair (just one in this case), find the _Uab_ matrix that describes the orientation of the internuclear vector.
 - Use Transpose(_Uab_) to build higher order rotation matrices, for the quadrupoles and higher.
 - Rotate the lab-frame multipoles to the quasi-internal frame, using these rotation matrices.
 - Form the intermediates for the _x_, _y_ and _z_ torques, from the quasi internal frame multipoles.
 - Construct the V^{ab}, V^{ba} and V^{ab,R} vectors, which are used for the energies, forces and torques.
 - Perform the appropriate dot products to form the energy, and the intermediates for the forces and torques.
 - Build the F vector from those intermediates, and rotate back to the lab frame, using the _Uab_ matrix.
 - Accumulate the torque on centers _a_ and _b_.
 - After the pair loop, map the torque on each atom back to forces on that atom and its anchors.
