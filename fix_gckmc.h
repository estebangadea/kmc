/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(gckmc,FixGCkMC)

#else

#ifndef LMP_FIX_GCKMC_H
#define LMP_FIX_GCKMC_H

#include <stdio.h>
#include "fix.h"
//#include "pointers.h"   // added by Jibao, according to Matias' 2012 version

namespace LAMMPS_NS {

class FixGCkMC : public Fix {
 public:
  FixGCkMC(class LAMMPS *, int, char **);
  ~FixGCkMC();
  int setmask();
  void init();
  //void pre_exchange();
  void END_OF_STEP();
  void attempt_atomic_freaction(int, double *);
  void attempt_atomic_breaction(int, double *);
  void extract_surface_position(double *, int);

  void update_gas_atoms_list();
  void update_reactive_atoms_list();
  void update_product_atoms_list();
  void update_surface_atoms_list();

  double compute_vector(int);
  double memory_usage();
  //void write_restart(FILE *);
  //void restart(char *);

    //void create_gaslist(); // from Matias' version; added by Jibao

 private:
    //int pressflag;      // 0=no 1=yes        // from Matias' version; added by Jibao

    FILE *fp;             // from Matias' version; added by Jibao

    //class PairHybrid *pairhybrid; // from Matias' version; added by Jibao



  int molecule_group,molecule_group_bit;
  int molecule_group_inversebit;
  int exclusion_group,exclusion_group_bit;
  int ngcmc,nevery,seed;
  int reactive_type, product_type, surf_type;
  int nreactions;
  int ngas;                 // # of gas atoms on all procs
  int nreact, nprod, nsurf;
  int ngas_local;           // # of gas atoms on this proc
  int nreact_local, nprod_local, nsurf_local; //Esteban
  int ngas_before;          // # of gas atoms on procs < this proc
  int nreact_before, nprod_before, nsurf_before;
  int mode;                 // ATOM or MOLECULE
  int regionflag;           // 0 = anywhere in box, 1 = specific region
  class Region *iregion;              // gcmc region
  char *idregion;           // gcmc region id
  int natoms_per_molecule;  // number of atoms in each gas molecule

  int groupbitall;          // group bitmask for inserted atoms

  double nfreaction_attempts;
  double nfreaction_successes;
  double nbreaction_attempts;
  double nbreaction_successes;
  int gcmc_nmax;
  double gas_mass;
  double reservoir_temperature;
  double tfac_insert;
  double displace;
  double beta,zz,sigma,volume;
  double kfreact, kbreact, potential, preexp, electrode_radi, electrode_h;//Esteban
  double center[3]; //Esteban
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double region_xlo,region_xhi,region_ylo,region_yhi,region_zlo,region_zhi;
  double region_volume;
  double *sublo,*subhi;
  int *local_gas_list;
  int *local_react_list;
  int *local_prod_list;
  int *local_surf_list;
  double **cutsq;
  double **atom_coord;
  double *surf_x;
  imageint imagezero;
  int imol,nmol;

  double energy_intra;
  class Molecule **onemols;
  class Pair *pair;

  class RanPark *random_equal;
  class RanPark *random_unequal;

  class Atom *model_atom;

  double **coords;
  imageint *imageflags;
  int triclinic;                         // 0 = orthog box, 1 = triclinic

  class Compute *c_pe;

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix gcmc does not (yet) work with atom_style template

Self-explanatory.

E: Fix gcmc region does not support a bounding box
Not all regions represent bounded volumes.  You cannot use
such a region with the fix gcmc command.

E: Fix gcmc region cannot be dynamic

Only static regions can be used with fix gcmc.

E: Fix gcmc region extends outside simulation box

Self-explanatory.

E: Fix gcmc molecule must have coordinates

The defined molecule does not specify coordinates.

E: Fix gcmc molecule must have atom types

The defined molecule does not specify atom types.

E: Atom type must be zero in fix gcmc mol command

Self-explanatory.

E: Fix gcmc molecule has charges, but atom style does not

Self-explanatory.

E: Fix gcmc molecule template ID must be same as atom_style template ID

When using atom_style template, you cannot insert molecules that are
not in that template.

E: Fix gcmc atom has charge, but atom style does not

Self-explanatory.

E: Cannot use fix gcmc shake and not molecule

Self-explanatory.

E: Molecule template ID for fix gcmc does not exist

Self-explanatory.

W: Molecule template for fix gcmc has multiple molecules

The fix gcmc command will only create molecules of a single type,
i.e. the first molecule in the template.

E: Region ID for fix gcmc does not exist

Self-explanatory.

W: Fix gcmc using full_energy option

Fix gcmc has automatically turned on the full_energy option since it
is required for systems like the one specified by the user. User input
included one or more of the following: kspace, triclinic, a hybrid
pair style, an eam pair style, or no "single" function for the pair
style.

E: Invalid atom type in fix gcmc command

The atom type specified in the gcmc command does not exist.

E: Fix gcmc cannot exchange individual atoms belonging to a molecule

This is an error since you should not delete only one atom of a
molecule.  The user has specified atomic (non-molecular) gas
exchanges, but an atom belonging to a molecule could be deleted.

E: All mol IDs should be set for fix gcmc group atoms

The molecule flag is on, yet not all molecule ids in the fix group
have been set to non-zero positive values by the user. This is an
error since all atoms in the fix gcmc group are eligible for deletion,
rotation, and translation and therefore must have valid molecule ids.

E: Fix gcmc molecule command requires that atoms have molecule attributes

Should not choose the gcmc molecule feature if no molecules are being
simulated. The general molecule flag is off, but gcmc's molecule flag
is on.

E: Fix gcmc shake fix does not exist

Self-explanatory.

E: Fix gcmc and fix shake not using same molecule template ID

Self-explanatory.

E: Cannot use fix gcmc in a 2d simulation

Fix gcmc is set up to run in 3d only. No 2d simulations with fix gcmc
are allowed.

E: Could not find fix gcmc exclusion group ID

Self-explanatory.

E: Could not find fix gcmc rotation group ID

Self-explanatory.

E: Illegal fix gcmc gas mass <= 0

The computed mass of the designated gas molecule or atom type was less
than or equal to zero.

E: Cannot do GCMC on atoms in atom_modify first group

This is a restriction due to the way atoms are organized in a list to
enable the atom_modify first command.

E: Could not find specified fix gcmc group ID

Self-explanatory.

E: Fix gcmc put atom outside box

This should not normally happen.  Contact the developers.

E: Fix gcmc ran out of available molecule IDs

See the setting for tagint in the src/lmptype.h file.

E: Fix gcmc ran out of available atom IDs

See the setting for tagint in the src/lmptype.h file.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

*/
