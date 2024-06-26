/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   On February 9, 2015, added 'CrackAniso' style for body centred cubic materials by Praveenkumar Hiremath (TU Bergakademie Freiberg | Computational Material Science (Years 2013-15) | www.tu-freiberg.de |)
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(displace_atoms,DisplaceAtoms)

#else
                       
#ifndef LMP_DISPLACE_ATOMS_H
#define LMP_DISPLACE_ATOMS_H

#include "pointers.h"

namespace LAMMPS_NS {

class DisplaceAtoms : protected Pointers {
 public:
  DisplaceAtoms(class LAMMPS *);
  void command(int, char **);
  void rotation(double o,double p,double q,double u,double v,double w,double s11,double s12,double s44,double s11_New,double s12_New,double s22_New,double s16_New,double s26_New,double s66_New,double A[]);
  void complex_add(double sum[],double Cmplx1[],double Cmplx2[]);
  void complex_subtract(double sub[],double Cmplx1[],double Cmplx2[]);
  void complex_multi(double mul[],double Cmplx1[],double Cmplx2[]);
  void complex_mul_inverse(double inv[],double Cmplx[]);
  void complex_square(double square[],double Cmplx1[]);
  void complex_root(double sq_root[],double Cmplx1[]);



 private:
  int scaleflag;
  void options(int, char **);
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Displace_atoms command before simulation box is defined

The displace_atoms command cannot be used before a read_data,
read_restart, or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot displace_atoms after reading restart file with per-atom info

This is because the restart file info cannot be migrated with the
atoms.  You can get around this by performing a 0-timestep run which
will assign the restart file info to actual atoms.

E: Could not find displace_atoms group ID

Group ID used in the displace_atoms command does not exist.

E: Invalid displace_atoms rotate axis for 2d

Axis must be in z direction.

E: Zero length rotation vector with displace_atoms

Self-explanatory.

W: Lost atoms via displace_atoms: original %ld current %ld

The command options you have used caused atoms to be lost.

*/
