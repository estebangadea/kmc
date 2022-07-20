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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier, Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_gckmc.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
#include "molecule.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "compute.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "thermo.h"
#include "output.h"
#include "neighbor.h"
#include <iostream>
#include "pair_hybrid_overlay.h"    // added by Jibao
#include "pair_sw.h"        // added by Jibao
//#include "pair_sw0.h"        // added by Jibao
#include "pair_hybrid.h"        // added by Jibao


using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{ATOM,MOLECULE};

/* ---------------------------------------------------------------------- */

FixGCkMC::FixGCkMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Esteban: Illegal fix gcmc command");

  if (atom->molecular == 2)
    error->all(FLERR,"Fix gcmc does not (yet) work with atom_style template");

  dynamic_group_allow = 1;

  vector_flag = 1;
  //size_vector = 8; // commented out by Jibao
    size_vector = 14;    // added by Jibao according to Matias' 2012 lammps version, to output energyout
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args


  reservoir_temperature = force->numeric(FLERR,arg[3]);
  reactive_type = force->inumeric(FLERR,arg[4]);
  product_type = force->inumeric(FLERR,arg[5]);
  preexp = force->numeric(FLERR,arg[6]);
  potential = force->numeric(FLERR,arg[7]);
  electrode_radi = force->numeric(FLERR,arg[8]);
  electrode_h = force->numeric(FLERR,arg[9]);
  seed = force->inumeric(FLERR,arg[10]);

 //Esteban: reactive_type, product_type, E, region, nreactions

  if (seed <= 0) error->all(FLERR,"Illegal fix gcmc command");
  if (reservoir_temperature < 0.0)
    error->all(FLERR,"Illegal fix gcmc command");

    //molflag = 0; // variable in 2012 verion // Jibao
    pairflag = 0; // added by Jibao. from Matias
    //pressflag=0;    // added by Jibao. from Matias
    regionflag=0;   // added by Jibao. from Matias

  // read options from end of input line

  options(narg-12,&arg[12]);

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp,seed);

  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi =
    region_zlo = region_zhi = 0.0;
  if (regionflag) {
    if (domain->regions[iregion]->bboxflag == 0)
      error->all(FLERR,"Fix gcmc region does not support a bounding box");
    if (domain->regions[iregion]->dynamic_check())
      error->all(FLERR,"Fix gcmc region cannot be dynamic");

    region_xlo = domain->regions[iregion]->extent_xlo;
    region_xhi = domain->regions[iregion]->extent_xhi;
    region_ylo = domain->regions[iregion]->extent_ylo;
    region_yhi = domain->regions[iregion]->extent_yhi;
    region_zlo = domain->regions[iregion]->extent_zlo;
    region_zhi = domain->regions[iregion]->extent_zhi;

    if (region_xlo < domain->boxlo[0] || region_xhi > domain->boxhi[0] ||
        region_ylo < domain->boxlo[1] || region_yhi > domain->boxhi[1] ||
        region_zlo < domain->boxlo[2] || region_zhi > domain->boxhi[2])
      error->all(FLERR,"Fix gcmc region extends outside simulation box");

    // estimate region volume using MC trials

    double coord[3];
    int inside = 0;
    int attempts = 10000000;
    for (int i = 0; i < attempts; i++) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      if (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) != 0)
        inside++;
    }

    double max_region_volume = (region_xhi - region_xlo)*
     (region_yhi - region_ylo)*(region_zhi - region_zlo);

    region_volume = max_region_volume*static_cast<double> (inside)/
     static_cast<double> (attempts);
  }

  // error check and further setup for mode = MOLECULE

  if (mode == MOLECULE) {
    if (onemols[imol]->xflag == 0)
      error->all(FLERR,"Fix gcmc molecule must have coordinates");
    if (onemols[imol]->typeflag == 0)
      error->all(FLERR,"Fix gcmc molecule must have atom types");
    if (ngcmc_type != 0)
      error->all(FLERR,"Atom type must be zero in fix gcmc mol command");
    if (onemols[imol]->qflag == 1 && atom->q == NULL)
      error->all(FLERR,"Fix gcmc molecule has charges, but atom style does not");

    if (atom->molecular == 2 && onemols != atom->avec->onemols)
      error->all(FLERR,"Fix gcmc molecule template ID must be same "
                 "as atom_style template ID");
    onemols[imol]->check_attributes(0);
  }

  if (charge_flag && atom->q == NULL)
    error->all(FLERR,"Fix gcmc atom has charge, but atom style does not");

  if (shakeflag && mode == ATOM)
    error->all(FLERR,"Cannot use fix gcmc shake and not molecule");

  // setup of coords and imageflags array

  if (mode == ATOM) natoms_per_molecule = 1;
  else natoms_per_molecule = onemols[imol]->natoms;
  memory->create(coords,natoms_per_molecule,3,"gcmc:coords");
  memory->create(imageflags,natoms_per_molecule,"gcmc:imageflags");
  memory->create(atom_coord,natoms_per_molecule,3,"gcmc:atom_coord");

  // compute the number of MC cycles that occur nevery timesteps

  //ncycles = nexchanges + nmcmoves + nreactions; //Esteban: Agregar nreactions

  // set up reneighboring
  force_reneighbor = 0;
  global_freq = nevery;
  // zero out counters

  ntranslation_attempts = 0.0;
  ntranslation_successes = 0.0;
  nrotation_attempts = 0.0;
  nrotation_successes = 0.0;
  ndeletion_attempts = 0.0;
  ndeletion_successes = 0.0;
  ninsertion_attempts = 0.0;
  ninsertion_successes = 0.0;
  nfreaction_attempts = 0.0;
  nfreaction_successes = 0.0;
  nbreaction_attempts = 0.0;
  nbreaction_successes = 0.0;


  //Esteban: nfreaction_attempts, nbreaction_attempts, nfreaction_successes, nbreaction_successes

    energyout=0.0;  // Matias

  gcmc_nmax = 0;
  local_gas_list = NULL;
  local_react_list = NULL; //Esteban
  local_prod_list = NULL;  //Esteban

   // if (comm->me == 0) printf("End of FixGCkMC::FixGCkMC()\n");
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixGCkMC::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix gcmc command");

  // defaults

  mode = ATOM;
  max_rotation_angle = 10*MY_PI/180;
  regionflag = 0;
  iregion = -1;
  region_volume = 0;
  max_region_attempts = 1000;
  molecule_group = 0;
  molecule_group_bit = 0;
  molecule_group_inversebit = 0;
  exclusion_group = 0;
  exclusion_group_bit = 0;
  pressure_flag = false;
  pressure = 0.0;
  fugacity_coeff = 1.0;
  shakeflag = 0;
  charge = 0.0;
  charge_flag = false;
  full_flag = false;
  idshake = NULL;
  ngroups = 0;
  int ngroupsmax = 0;
  groupstrings = NULL;
  ngrouptypes = 0;
  int ngrouptypesmax = 0;
  grouptypestrings = NULL;
  grouptypes = NULL;
  grouptypebits = NULL;
  energy_intra = 0.0;
  tfac_insert = 1.0;

  int iarg = 0;
  while (iarg < narg) {
  if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      imol = atom->find_molecule(arg[iarg+1]);
      if (imol == -1)
        error->all(FLERR,"Molecule template ID for fix gcmc does not exist");
      if (atom->molecules[imol]->nset > 1 && comm->me == 0)
        error->warning(FLERR,"Molecule template for "
                       "fix gcmc has multiple molecules");
      mode = MOLECULE;
      onemols = atom->molecules;
      nmol = onemols[imol]->nset;
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix gcmc does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxangle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      max_rotation_angle = force->numeric(FLERR,arg[iarg+1]);
      max_rotation_angle *= MY_PI/180;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pair") == 0) { // added by Jibao. from Matias
        if (iarg+2 > narg) error->all(FLERR,"Illegal fix GCMC command");
        if (strcmp(arg[iarg+1],"lj/cut") == 0) pairflag = 0;
        else if (strcmp(arg[iarg+1],"Stw") == 0) pairflag = 1;
        else error->all(FLERR,"Illegal fix evaporate command");
        iarg += 2;
    }   // Matias
    else if (strcmp(arg[iarg],"pressure") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      pressure = force->numeric(FLERR,arg[iarg+1]);
        pressure = pressure * 100.0;    // added by Jibao, according to Matias' code
      pressure_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"fugacity_coeff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      fugacity_coeff = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      charge = force->numeric(FLERR,arg[iarg+1]);
      charge_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"shake") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      int n = strlen(arg[iarg+1]) + 1;
      delete [] idshake;
      idshake = new char[n];
      strcpy(idshake,arg[iarg+1]);
      shakeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"full_energy") == 0) {
      full_flag = true;
      iarg += 1;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      if (ngroups >= ngroupsmax) {
	ngroupsmax = ngroups+1;
	groupstrings = (char **)
	  memory->srealloc(groupstrings,
			   ngroupsmax*sizeof(char *),
			   "fix_gcmc:groupstrings");
      }
      int n = strlen(arg[iarg+1]) + 1;
      groupstrings[ngroups] = new char[n];
      strcpy(groupstrings[ngroups],arg[iarg+1]);
      ngroups++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"grouptype") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix gcmc command");
      if (ngrouptypes >= ngrouptypesmax) {
	ngrouptypesmax = ngrouptypes+1;
	grouptypes = (int*) memory->srealloc(grouptypes,ngrouptypesmax*sizeof(int),
			 "fix_gcmc:grouptypes");
	grouptypestrings = (char**)
	  memory->srealloc(grouptypestrings,
			   ngrouptypesmax*sizeof(char *),
			   "fix_gcmc:grouptypestrings");
      }
      grouptypes[ngrouptypes] = atoi(arg[iarg+1]);
      int n = strlen(arg[iarg+2]) + 1;
      grouptypestrings[ngrouptypes] = new char[n];
      strcpy(grouptypestrings[ngrouptypes],arg[iarg+2]);
      ngrouptypes++;
      iarg += 3;
    } else if (strcmp(arg[iarg],"intra_energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      energy_intra = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tfac_insert") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix gcmc command");
      tfac_insert = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix gcmc command");
  }
     // if (comm->me == 0) printf("End of FixGCkMC::options()\n");
}

/* ---------------------------------------------------------------------- */

FixGCkMC::~FixGCkMC()
{
   // printf("FixGCkMC()");
  if (regionflag) delete [] idregion;
  delete random_equal;
  delete random_unequal;

    //delete region_insert;   // from Matias; deleted by Jibao

  memory->destroy(local_gas_list);
  memory->destroy(local_react_list);
  memory->destroy(local_prod_list);
  memory->destroy(atom_coord);
  memory->destroy(coords);
  memory->destroy(imageflags);

  delete [] idshake;

  if (ngroups > 0) {
    for (int igroup = 0; igroup < ngroups; igroup++)
      delete [] groupstrings[igroup];
    memory->sfree(groupstrings);
  }

  if (ngrouptypes > 0) {
    memory->destroy(grouptypes);
    memory->destroy(grouptypebits);
    for (int igroup = 0; igroup < ngrouptypes; igroup++)
      delete [] grouptypestrings[igroup];
    memory->sfree(grouptypestrings);
  }
   // if (comm->me == 0) printf("End of FixGCkMC::~FixGCkMC()\n");
}

/* ---------------------------------------------------------------------- */

int FixGCkMC::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGCkMC::init()
{
   // if (comm->me == 0) printf("Begins: FixGCkMC::init()\n");

  triclinic = domain->triclinic;

  // decide whether to switch to the full_energy option

  if (!full_flag) {
    if ((force->kspace) ||
        (force->pair == NULL) ||
        (force->pair->single_enable == 0) ||
        (force->pair_match("hybrid",0)) ||
        (force->pair_match("eam",0))
	) {
      full_flag = true;

        //if (comm->me == 0) printf("Begins: inside if (!full_flag){}: FixGCkMC::init()\n");

        //if (comm->me == 0) printf("pairflag = %d\n",pairflag);

        if (pairflag) { // added by Jibao
            full_flag = false;  // added by Jibao
        }   // added by Jibao

      if (comm->me == 0 && full_flag == true) // modified by Jibao
          error->warning(FLERR,"Fix gcmc using full_energy option");
    }
  }

    //if (comm->me == 0) printf("Begins 2: FixGCkMC::init()\n");

  if (full_flag) {
    char *id_pe = (char *) "thermo_pe";
    int ipe = modify->find_compute(id_pe);
    c_pe = modify->compute[ipe];
  }

  int *type = atom->type;

  if (mode == ATOM) {
    if (product_type <= 0 || product_type > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix gcmc command");
    if (reactive_type <= 0 || reactive_type > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix gcmc command");
  }

    //if (comm->me == 0) printf("Begins 3: FixGCkMC::init()\n");

  // if mode == ATOM, warn if any deletable atom has a mol ID

  if ((mode == ATOM) && atom->molecule_flag) {
      /*
      if (comm->me == 0) {
          printf("Inside if ((mode == ATOM)): FixGCkMC::init()\n");
          printf("atom->molecule_flag = %d\n",atom->molecule_flag);
      }
      */
    tagint *molecule = atom->molecule;
    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if (type[i] == reactive_type)
        if (molecule[i]) flag = 1;

      //if (comm->me == 0) printf("Inside if ((mode == ATOM)) 2: FixGCkMC::init()\n");

    int flagall;

      //printf("comm->me = %d, flag = %d, flagall = %d, before MPI_ALLreduce()\n",comm->me,flag,flagall);

      //error->all(FLERR,"Kao 0 !!!!");    // added by Jibao

    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

      //error->all(FLERR,"Kao 1 !!!!");    // added by Jibao

      //if (comm->me == 0) printf("Inside if ((mode == ATOM)) 3: FixGCkMC::init()\n");
      //if (comm->me == 0) printf("flag = %d, flagall = %d, after MPI_ALLreduce()\n",flag,flagall);

      if (flagall && comm->me == 0) {
          //if (comm->me == 0) printf("Inside if if (flagall && comm->me == 0): FixGCkMC::init()\n");    // added by Jibao
          //error->all(FLERR,"Kao 2 !!!!");    // added by Jibao
          error->all(FLERR,"Fix gcmc cannot exchange individual atoms belonging to a molecule");
      }
  }

    //if (comm->me == 0) printf("Begins 4: FixGCkMC::init()\n");

  // if mode == MOLECULE, check for unset mol IDs

  if (mode == MOLECULE) {
    tagint *molecule = atom->molecule;
    int *mask = atom->mask;
    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] == groupbit)
        if (molecule[i] == 0) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall && comm->me == 0)
      error->all(FLERR,
       "All mol IDs should be set for fix gcmc group atoms");
  }

    //if (comm->me == 0) printf("Begins 5: FixGCkMC::init()\n");

  if (((mode == MOLECULE) && (atom->molecule_flag == 0)) ||
      ((mode == MOLECULE) && (!atom->tag_enable || !atom->map_style)))
    error->all(FLERR,
               "Fix gcmc molecule command requires that "
               "atoms have molecule attributes");

  // if shakeflag defined, check for SHAKE fix
  // its molecule template must be same as this one

    //if (comm->me == 0) printf("Begins 6: FixGCkMC::init()\n");

  fixshake = NULL;
  if (shakeflag) {
    int ifix = modify->find_fix(idshake);
    if (ifix < 0) error->all(FLERR,"Fix gcmc shake fix does not exist");
    fixshake = modify->fix[ifix];
    int tmp;
    if (onemols != (Molecule **) fixshake->extract("onemol",tmp))
      error->all(FLERR,"Fix gcmc and fix shake not using "
                 "same molecule template ID");
  }

    //if (comm->me == 0) printf("Begins 7: FixGCkMC::init()\n");

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix gcmc in a 2d simulation");

  // create a new group for interaction exclusions

    //if (comm->me == 0) printf("Before 'create a new group for interaction exclusions': FixGCkMC::init()\n");

  if (full_flag || pairflag) {  // modified by Jibao; added "|| pairflag"
    char **group_arg = new char*[4];
    // create unique group name for atoms to be excluded
    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    sprintf(group_arg[0],"FixGCMC:gcmc_exclusion_group:%s",id);
    group_arg[1] = (char *) "subtract";
    group_arg[2] = (char *) "all";
    group_arg[3] = (char *) "all";
    group->assign(4,group_arg);
    exclusion_group = group->find(group_arg[0]);
    if (exclusion_group == -1)
      error->all(FLERR,"Could not find fix gcmc exclusion group ID");
    exclusion_group_bit = group->bitmask[exclusion_group];

    // neighbor list exclusion setup
    // turn off interactions between group all and the exclusion group

    int narg = 4;
    char **arg = new char*[narg];;
    arg[0] = (char *) "exclude";
    arg[1] = (char *) "group";
    arg[2] = group_arg[0];
    arg[3] = (char *) "all";
    neighbor->modify_params(narg,arg);
    delete [] group_arg[0];
    delete [] group_arg;
    delete [] arg;
  }

  // create a new group for temporary use with selected molecules

  if (mode == MOLECULE) {
    char **group_arg = new char*[3];
    // create unique group name for atoms to be rotated
    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    sprintf(group_arg[0],"FixGCMC:rotation_gas_atoms:%s",id);
    group_arg[1] = (char *) "molecule";
    char digits[12];
    sprintf(digits,"%d",-1);
    group_arg[2] = digits;
    group->assign(3,group_arg);
    molecule_group = group->find(group_arg[0]);
    if (molecule_group == -1)
      error->all(FLERR,"Could not find fix gcmc rotation group ID");
    molecule_group_bit = group->bitmask[molecule_group];
    molecule_group_inversebit = molecule_group_bit ^ ~0;
    delete [] group_arg[0];
    delete [] group_arg;
  }

  // get all of the needed molecule data if mode == MOLECULE,
  // otherwise just get the gas mass

  if (mode == MOLECULE) {

    onemols[imol]->compute_mass();
    onemols[imol]->compute_com();
    gas_mass = onemols[imol]->masstotal;
    for (int i = 0; i < onemols[imol]->natoms; i++) {
      onemols[imol]->x[i][0] -= onemols[imol]->com[0];
      onemols[imol]->x[i][1] -= onemols[imol]->com[1];
      onemols[imol]->x[i][2] -= onemols[imol]->com[2];
    }

  } else //gas_mass = atom->mass[ngcmc_type];

  //if (gas_mass <= 0.0)
  //  error->all(FLERR,"Illegal fix gcmc gas mass <= 0");

  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot do GCMC on atoms in atom_modify first group");
  }

    //if (comm->me == 0) printf("Before 'compute beta, lambda, sigma, and the zz factor': FixGCkMC::init()\n");

  // compute beta, lambda, sigma, and the zz factor

  beta = 1.0/(1.38065e-23*reservoir_temperature);

  kfreact = preexp*exp(0.5*1*1.6022e-19*beta*potential); //Esteban: Potencial estandard de la oxidacion del agua en V.
  kbreact = preexp*exp(-0.5*1*1.6022e-19*beta*potential);
  center[0] = (domain->boxhi[0]-domain->boxlo[0])/2;
  center[1] = (domain->boxhi[1]-domain->boxlo[1])/2;
  center[2] = (domain->boxhi[2]-domain->boxlo[2])/2;

//  double lambda = sqrt(force->hplanck*force->hplanck/
//                       (2.0*MY_PI*gas_mass*force->mvv2e*
//                        force->boltz*reservoir_temperature));
//  sigma = sqrt(force->boltz*reservoir_temperature*tfac_insert/gas_mass/force->mvv2e);

//    if (!pressure_flag) // using pressure_flag to replace pressflag defined by Matias
//        zz = exp(beta*chemical_potential)/(pow(lambda,3.0));
//    else if(pressure_flag) {
//        zz = pressure/(force->boltz*4.184*reservoir_temperature*1000/6.02e23)/(1e30);   // from Matias; need to check the meaning
        //zz = pressure*fugacity_coeff*beta/force->nktv2p;    // from the original expression of 2015 version of lammps
//    }


  //zz = exp(beta*chemical_potential)/(pow(lambda,3.0)); // commented out by Jibao
  //if (pressure_flag) zz = pressure*fugacity_coeff*beta/force->nktv2p; // commented out by Jibao

    if (comm->me==0) { // added by Jibao; from Matias' version of lammps
        //printf("zz factor equals %e\n",zz); // added by Jibao; from Matias' version of lammps
        printf("regionflag equals %i\n",regionflag); // added by Jibao; from Matias' version of lammps
        printf("pressure_flag equals %i\n",pressure_flag); // added by Jibao; from Matias' version of lammps
        printf("pairflag equals %i\n",pairflag); // added by Jibao; from Matias' version of lammps
    } // added by Jibao; from Matias' version of lammps

  imagezero = ((imageint) IMGMAX << IMG2BITS) |
             ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // construct group bitmask for all new atoms
  // aggregated over all group keywords

  groupbitall = 1 | groupbit;
  for (int igroup = 0; igroup < ngroups; igroup++) {
    int jgroup = group->find(groupstrings[igroup]);
    if (jgroup == -1)
      error->all(FLERR,"Could not find specified fix gcmc group ID");
    groupbitall |= group->bitmask[jgroup];
  }

  // construct group type bitmasks
  // not aggregated over all group keywords

  if (ngrouptypes > 0) {
    memory->create(grouptypebits,ngrouptypes,"fix_gcmc:grouptypebits");
    for (int igroup = 0; igroup < ngrouptypes; igroup++) {
      int jgroup = group->find(grouptypestrings[igroup]);
      if (jgroup == -1)
	error->all(FLERR,"Could not find specified fix gcmc group ID");
      grouptypebits[igroup] = group->bitmask[jgroup];
    }
  }

    //printf("End of FixGCkMC::init()\n");

}

/* ----------------------------------------------------------------------
   attempt Monte Carlo translations, rotations, insertions, and deletions
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixGCkMC::pre_exchange()
{
  // just return if should not be called on this timestep
 //if (comm->me == 0) printf("Begin of FixGCkMC::pre_exchange()\n");
  //if (next_reneighbor != update->ntimestep) return;

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  if (triclinic) {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } else {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }

  if (regionflag) volume = region_volume;
  else volume = domain->xprd * domain->yprd * domain->zprd;

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  update_gas_atoms_list();
  update_product_atoms_list();
  update_reactive_atoms_list();

  //printf("nprod=%d nreact=%d\n", nprod_local, nreact_local);

  if (full_flag) {
    error->all(FLERR,"gkcmc does not allow full energy");

  } else {

    if (mode == MOLECULE) {
      error->all(FLERR,"gkcmc does not allow mode MOLECULE");
    } else {
        for (int i = 0; i < nreact; i++){
            attempt_atomic_freaction(i);
        }
        for (int i = 0; i < nprod; i++){
            attempt_atomic_breaction(i);
         }
     }
            //Esteban: Agregar la probabilidad de una freaction o breaction

      //domain->pbc();    // added by Jibao; to prevent the error: ERROR on proc 0: Bond atoms 4205 4209 missing on proc 0 at step 285 (../neigh_bond.cpp:196)
      //comm->exchange(); // added by Jibao; to prevent the error: ERROR on proc 0: Bond atoms 4205 4209 missing on proc 0 at step 285 (../neigh_bond.cpp:196)
  }
 //if (comm->me == 0) printf("End of FixGCkMC::pre_exchange()\n");
}

void FixGCkMC::attempt_atomic_freaction(int i)
{
  nfreaction_attempts += 1.0;

  int success = 0;

  if ((i >= nreact_before) &&
      (i < nreact_before + nreact_local)){
  int ilocal= i - nreact_before;
  //printf("freaction try\n");
  int j = local_react_list[ilocal];
  int tstep = update->dt;

    double **x = atom->x;
    double kvel;
    int *type = atom->type;
    //double energy_before = energy(i, reactive_type,-1,x[i]);
    //double energy_after = energy(i, product_type,-1,x[i]);
    double r_xy = (x[j][0]-center[0])*(x[j][0]-center[0])
        +(x[j][1]-center[1])*(x[j][1]-center[1]);

    if (r_xy < electrode_radi*electrode_radi){
        kvel = exp(electrode_h-x[j][2])*kfreact;
    }
    else{
        double r_xyz = sqrt((sqrt(r_xy)-electrode_radi)*(sqrt(r_xy)-electrode_radi)+
            (x[j][2]-electrode_h)*(x[j][2]-electrode_h));
        kvel = exp(-r_xyz)*kfreact;
    }

    if (random_unequal->uniform() <
        1-exp(-kvel*tstep)) {
            type[j] = product_type;
            success = 1;
            printf("freaction: x=%f  y=%f  p=%e\n", x[j][0],x[j][1],1-exp(-kvel*5));
    }
}


    int success_all = 0;
    //printf("antes del MPI\n");
    MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);
    //printf("despues del MPI\n");
    if (success_all) {
        update_gas_atoms_list();
        update_reactive_atoms_list();
        update_product_atoms_list();
        nfreaction_successes += 1;

        atom->nghost = 0;
        comm->borders();

    }
    //printf("freaction end");
}
 /* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCkMC::attempt_atomic_breaction(int i)
{
  nbreaction_attempts += 1.0;

  int success = 0;

  if ((i >= nprod_before) &&
      (i < nprod_before + nprod_local)){

  //printf("breaction try\n");
    int ilocal = i - nprod_before;
    int j = local_prod_list[ilocal];
    int tstep = update->dt;

    double **x = atom->x;
    double kvel;
    int *type = atom->type;
    //double energy_before = energy(i, reactive_type,-1,x[i]);
    //double energy_after = energy(i, product_type,-1,x[i]);
    double r_xy = (x[j][0]-center[0])*(x[j][0]-center[0])
        +(x[j][1]-center[1])*(x[j][1]-center[1]);

    if (r_xy < electrode_radi*electrode_radi){
        kvel = exp(electrode_h-x[j][2])*kbreact;
    }
    else{
        double r_xyz = sqrt((sqrt(r_xy)-electrode_radi)*(sqrt(r_xy)-electrode_radi)+
            (x[j][2]-electrode_h)*(x[j][2]-electrode_h));
        kvel = exp(-r_xyz)*kbreact;
    }

    if (random_unequal->uniform() <
        1-exp(-kvel*tstep)) {
            type[j] = reactive_type;
            success = 1;
            printf("breaction: x=%f  y=%f  p=%e\n", x[j][0],x[j][1],1-exp(-kvel*5));
  }
}

    int success_all = 0;
    MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);

    if (success_all) {
        update_gas_atoms_list();
        update_reactive_atoms_list();
        update_product_atoms_list();
        nbreaction_successes += 1;

        atom->nghost = 0;
        comm->borders();
    }
}


/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double FixGCkMC::energy(int i, int itype, tagint imolecule, double *coord)
{
    //printf("comm->me = %d: Beginning of FixGCkMC::energy()\n",comm->me);
    double delx,dely,delz,rsq;
    double rsq1, rsq2;  // added by Jibao; for Stw_GCMC
    int ietype,jetype,ketype,ijparam,ikparam,ijkparam;  // added by Jibao; ietype: element type of itype
    int jtype,ktype;
    double delr1[3],delr2[3],fj[3],fk[3];   // added by Jibao; for Stw_GCMC
    int jj,kk;  // added by Jibao; for Stw_GCMC

    double **x = atom->x;
    int *type = atom->type;
    tagint *molecule = atom->molecule;
    int nall = atom->nlocal + atom->nghost;

    //printf("comm->me = %d: nall= %d, atom->nlocal= %d, atom->nghost= %d\n",comm->me,nall,atom->nlocal,atom->nghost);

    pair = force->pair;
    cutsq = force->pair->cutsq;

    double fpair = 0.0;
    double factor_coul = 1.0;
    double factor_lj = 1.0;

    double total_energy = 0.0;
    double twobodyeng = 0.0;    // added by Jibao; for Stw_GCMC
    double tmp2body = 0.0;      // added by Jibao; for Stw_GCMC
    double threebodyeng = 0.0;  // added by Jibao; for Stw_GCMC
    double tmp3body = 0.0;      // added by Jibao; for Stw_GCMC

    char *pair_style;   // added by Jibao
    pair_style = force->pair_style; // added by Jibao
    int eflag = 1;  // added by Jibao; for twobody() and threebody()

    //printf("comm->me = %d: 2 of FixGCkMC::energy()\n",comm->me);

    if (!strcmp(pair_style,"hybrid")) {
        int ***map_substyle = (int ***) pair->returnmap_substyle();   // added by Jibao; return list of sub-styles itype,jtype points to

        Pair **styles = (Pair **) pair->returnstyles();     // added by Jibao; return list of Pair style classes

        char **keywords = (char **) pair->returnkeywords();   // added by Jibao; return style name of each Pair style

        //printf("comm->me = %d: 3 of FixGCkMC::energy()\n",comm->me);

        for(int j=0; j < nall; j++){

            //printf("comm->me = %d: 3.1 of FixGCkMC::energy()\n",comm->me);

            if (i == j) continue;
            if (mode == MOLECULE)
                if (imolecule == molecule[j]) continue;

            //printf("comm->me = %d: 3.2 of FixGCkMC::energy()\n",comm->me);

            delx = coord[0] - x[j][0];
            dely = coord[1] - x[j][1];
            delz = coord[2] - x[j][2];
            rsq = delx*delx + dely*dely + delz*delz;

            //printf("comm->me = %d: 3.3 of FixGCkMC::energy()\n",comm->me);

            //printf("comm->me = %d: j = %d\n",comm->me,j);

            jtype = type[j];

            //printf("comm->me = %d: type[j]= type[%d]= %d\n",comm->me,j,type[j]);

            //printf("comm->me = %d: 3.4 of FixGCkMC::energy()\n",comm->me);

            int **nmap = (int **) pair->returnnmap();
            //printf("comm->me = %d: nmap[%d][%d]= %d\n",comm->me,itype,jtype,nmap[itype][jtype]);

            if (nmap[itype][jtype] > 0) {
                //printf("comm->me = %d: itype= %d,jtype= %d,map_substyle[%d][%d][0]= %d\n",comm->me,itype,jtype,itype,jtype,map_substyle[itype][jtype][0]);

                //printf("comm->me = %d: keywords[map_substyle[%d][%d][0]] = %s\n",comm->me,itype,jtype,keywords[map_substyle[itype][jtype][0]]);

                int substyle = map_substyle[itype][jtype][0];

                if (strstr(keywords[substyle],"sw")) {

                    //if (!strcmp(keywords[substyle],"sw") || !strcmp(keywords[substyle],"sw0") || !strcmp(keywords[substyle],"sw/omp") || !strcmp(keywords[substyle],"sw0/omp")) {


                    //pair = force->pair_match("hybrid",0);

                    //printf("comm->me = %d: 4 of FixGCkMC::energy()\n",comm->me);

                    int *map = (int *) styles[substyle]->returnmap();
                    //int *map = (int *) pair->returnmap();

                    //printf("comm->me = %d: 4.1 of FixGCkMC::energy()\n",comm->me);
                    //printf("itype= %d in FixGCkMC::energy()\n",itype);

                    //for (int kao = 1; kao<=7; kao++) printf("map[%d] = %d\n",kao,map[kao]);


                    //printf("map[itype] = map[%d] = %d in FixGCkMC::energy()\n",itype,map[itype]);

                    ietype=map[itype];

                    //printf("comm->me = %d: 4.2 of FixGCkMC::energy()\n",comm->me);

                    LAMMPS_NS::Pair::Param *params = (LAMMPS_NS::Pair::Param *) styles[substyle]->returnparams();   // parameter set for an I-J-K interaction
                    //printf("comm->me = %d: 4.3 of FixGCkMC::energy()\n",comm->me);

                    int ***elem2param = (int ***) styles[substyle]->returnelem2param();
                    //printf("comm->me = %d: 4.4 of FixGCkMC::energy()\n",comm->me);

                    jetype=map[jtype];
                    //printf("comm->me = %d: 4.5 of FixGCkMC::energy()\n",comm->me);
                    //printf("comm->me = %d: itype= %d, ietype= %d, jtype= %d, map[jtype]= %d\n",comm->me,itype,ietype, jtype, map[jtype]);

                    ijparam = elem2param[ietype][jetype][jetype];
                    //printf("comm->me = %d: 4.6 of FixGCkMC::energy()\n",comm->me);

                    //printf("comm->me = %d: 5 of FixGCkMC::energy()\n",comm->me);

                    if (rsq < params[ijparam].cutsq){
                        styles[substyle]->twobody(&params[ijparam],rsq,fpair,eflag,tmp2body);
                        //printf("comm->me = %d: tmp2body = %f\n",comm->me,tmp2body);
                        total_energy+=tmp2body;
                        //printf("# %d %d %d %d %d %f %e\n",ietype,type[j],ietype,jetype,ijparam,params[ijparam].epsilon,tmp2body);
                    }
                    //printf("comm->me = %d: 6 of FixGCkMC::energy()\n",comm->me);
                } else {
                    //printf("comm->me = %d: 7 of FixGCkMC::energy()\n",comm->me);
                    if (rsq < cutsq[itype][jtype])
                        total_energy +=
                        styles[substyle]->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
                }
            }
            //printf("comm->me = %d: 8 of FixGCkMC::energy()\n",comm->me);
        }

        //printf("comm->me = %d: after twobody loop\n",comm->me);

        // only for sw, sw0 potential
        // first possibility ii!=i j==i k!=i or ii!=i; j!=i ; k==i
        for(int ii = 0;ii < nall;ii++){
            if (ii==i) continue;
            if (mode == MOLECULE)
                if (imolecule == molecule[ii]) continue;

            int **nmap = (int **) pair->returnnmap();
            //printf("first possibility: nmap[%d][%d]= %d in for(int ii = 0;ii < nall;ii++){}\n",jtype,itype,nmap[jtype][itype]);
            jtype = type[ii];

            if (nmap[jtype][itype] > 0) {

                int substyle_ji = map_substyle[jtype][itype][0];

                //printf("map_substyle[%d][%d][0]= %d in for(int ii = 0;ii < nall;ii++){}\n",jtype,itype,map_substyle[jtype][itype][0]);

                //if (comm->me == 0) printf("map_substyle[%d][%d][0] = %d, map_substyle[%d][%d][0] = %d\n",itype,jtype,map_substyle[itype][jtype][0],jtype,itype,map_substyle[jtype][itype][0]);

                if (!strcmp(keywords[substyle_ji],"sw")) {
                    int *map = (int *) styles[substyle_ji]->returnmap();
                    ietype= map[type[ii]];

                    LAMMPS_NS::Pair::Param *params = (LAMMPS_NS::Pair::Param *) styles[substyle_ji]->returnparams();   // parameter set for an I-J-K interaction
                    int ***elem2param = (int ***) styles[substyle_ji]->returnelem2param();

                    jetype=map[itype];

                    ijparam = elem2param[ietype][jetype][jetype];
                    delr1[0] = coord[0] - x[ii][0];
                    delr1[1] = coord[1] - x[ii][1];
                    delr1[2] = coord[2] - x[ii][2];
                    rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
                    if (rsq1 > params[ijparam].cutsq) continue;

                    for (int k = 0; k < nall; k++){
                        if (ii==k || k==i) continue;
                        if (mode == MOLECULE)
                            if (imolecule == molecule[ii]) continue;

                        //printf("first possibility: nmap[%d][%d]= %d in for(int k = 0; k < nall; k++){}\n",jtype,ktype,nmap[jtype][ktype]);
                        ktype = type[k];
                        if (nmap[jtype][ktype] > 0) {

                            int substyle_jk = map_substyle[jtype][ktype][0];

                            if (!strcmp(keywords[substyle_jk],"sw")) {  // if yes, substyle_jk == sybstyle_ji
                                ketype = map[type[k]];

                                delr2[0] = x[k][0] - x[ii][0];
                                delr2[1] = x[k][1] - x[ii][1];
                                delr2[2] = x[k][2] - x[ii][2];

                                ikparam = elem2param[ietype][ketype][ketype];
                                ijkparam = elem2param[ietype][jetype][ketype];

                                rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
                                if (rsq2 > params[ikparam].cutsq) continue;
                                styles[substyle_ji]->threebody(&params[ijparam],&params[ikparam],&params[ijkparam],rsq1,rsq2,delr1,delr2,fj,fk,eflag,tmp3body);

                                total_energy+=tmp3body;
                            }
                        }
                    }
                }
            }
        }

        //printf("comm->me = %d: after first possibility loop\n",comm->me);

        //second possibility ii==i j!=i k!=i

        for (int j = 0; j < nall; j++){
            if (i==j) continue;
            if (mode == MOLECULE)
                if (imolecule == molecule[j]) continue;

            int **nmap = (int **) pair->returnnmap();

            //printf("second possibility: nmap[%d][%d]= %d in for(int k = 0; k < nall; k++){}\n",itype,jtype,nmap[itype][jtype]);
            jtype = type[j];
            if (nmap[itype][jtype] > 0) {

                int substyle_ij = map_substyle[itype][jtype][0];

                if (!strcmp(keywords[substyle_ij],"sw")) {
                    int *map = (int *) styles[substyle_ij]->returnmap();
                    jetype= map[type[j]];
                    ietype=map[itype];

                    LAMMPS_NS::Pair::Param *params = (LAMMPS_NS::Pair::Param *) styles[substyle_ij]->returnparams();   // parameter set for an I-J-K interaction
                    int ***elem2param = (int ***) styles[substyle_ij]->returnelem2param();

                    ijparam = elem2param[ietype][jetype][jetype];
                    delr1[0] = x[j][0] - coord[0];
                    delr1[1] = x[j][1] - coord[1];
                    delr1[2] = x[j][2] - coord[2];
                    rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

                    if (rsq1 > params[ijparam].cutsq) continue;
                    for (int k = j; k < nall; k++) {
                        if (i==k || k == j ) continue;

                        ktype = type[k];

                        if (nmap[itype][ktype] > 0) {

                            int substyle_ik = map_substyle[itype][ktype][0];

                            if (!strcmp(keywords[substyle_ik],"sw")) {
                                ketype = map[type[k]];

                                ikparam = elem2param[ietype][ketype][ketype];
                                ijkparam = elem2param[ietype][jetype][ketype];
                                delr2[0] = x[k][0] - coord[0];
                                delr2[1] = x[k][1] - coord[1];
                                delr2[2] = x[k][2] - coord[2];
                                rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
                                if (rsq2 > params[ikparam].cutsq) continue;
                                styles[substyle_ij]->threebody(&params[ijparam],&params[ikparam],&params[ijkparam],rsq1,rsq2,delr1,delr2,fj,fk,eflag,tmp3body);
                                total_energy+=tmp3body;
                            }
                        }
                    }
                }
            }
        }
        //printf("comm->me = %d: after second possibility loop\n",comm->me);
        ////////// TOTAL ENERGY/////
        //total_energy=threebodyeng+twobodyeng;
    } else if (strstr(pair_style,"sw")) {
        //if (comm->me == 0) printf("energy() for Stw_GCMC part\n");
        int *map = (int *) pair->returnmap();
        ietype=map[itype];

        LAMMPS_NS::Pair::Param *params = (LAMMPS_NS::Pair::Param *) pair->returnparams();   // parameter set for an I-J-K interaction
        int ***elem2param = (int ***) pair->returnelem2param();


        for(int j = 0; j < nall; j++){
            //if (comm->me == 0) printf("j = %d in two-body part in Stw_GCMC()\n",j);    // added by Jibao
            if (i == j) continue;
            if (mode == MOLECULE)
                if (imolecule == molecule[j]) continue;
            //if (comm->me == 0) printf("x[%d][0] = %f,coord = %d\n",j,x[j][0],coord);// added by Jibao
            //if (comm->me == 0) printf("coord[0] = %f\n",coord[0]);// added by Jibao
            delx = coord[0] - x[j][0];
            dely = coord[1] - x[j][1];
            delz = coord[2] - x[j][2];
            //if (comm->me == 0) printf("j = %d in two-body part in Stw_GCMC()\n",j);    // added by Jibao
            rsq = delx*delx + dely*dely + delz*delz;
            jetype=map[type[j]];
            ijparam = elem2param[ietype][jetype][jetype];

            if (rsq < params[ijparam].cutsq){
                pair->twobody(&params[ijparam],rsq,fpair,eflag,tmp2body);        twobodyeng+=tmp2body;
                //printf("# %d %d %d %d %d %f %e\n",ietype,type[j],ietype,jetype,ijparam,params[ijparam].epsilon,tmp2body);
            }
        }

        //first possibility ii!=i j==i k!=i or ii!=i; j!=i ; k==i
        jetype=map[itype];
        for(int ii = 0;ii < nall;ii++){
            if (ii==i) continue;
            if (mode == MOLECULE)
                if (imolecule == molecule[ii]) continue;
            ietype= map[type[ii]];
            ijparam = elem2param[ietype][jetype][jetype];
            delr1[0] = coord[0] - x[ii][0];
            delr1[1] = coord[1] - x[ii][1];
            delr1[2] = coord[2] - x[ii][2];
            rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
            if (rsq1 > params[ijparam].cutsq) continue;
            for (int k = 0; k < nall; k++){
                if (ii==k || k==i) continue;
                ketype = map[type[k]];
                delr2[0] = x[k][0] - x[ii][0];
                delr2[1] = x[k][1] - x[ii][1];
                delr2[2] = x[k][2] - x[ii][2];
                ikparam = elem2param[ietype][ketype][ketype];
                ijkparam = elem2param[ietype][jetype][ketype];
                rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
                if (rsq2 > params[ikparam].cutsq) continue;
                pair->threebody(&params[ijparam],&params[ikparam],&params[ijkparam],rsq1,rsq2,delr1,delr2,fj,fk,eflag,tmp3body);
                threebodyeng+=tmp3body;
            }
        }
        //second possibility ii==i j!=i k!=i
        ietype=map[itype];
        for (int j = 0; j < nall; j++){
            if (i==j) continue;
            if (mode == MOLECULE)
                if (imolecule == molecule[j]) continue;
            jetype = map[type[j]];
            ijparam = elem2param[ietype][jetype][jetype];
            delr1[0] = x[j][0] - coord[0];
            delr1[1] = x[j][1] - coord[1];
            delr1[2] = x[j][2] - coord[2];
            rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
            if (rsq1 > params[ijparam].cutsq) continue;
            for (int k = j; k < nall; k++) {
                if (i==k || k == j ) continue;
                ketype = map[type[k]];
                ikparam = elem2param[ietype][ketype][ketype];
                ijkparam = elem2param[ietype][jetype][ketype];
                delr2[0] = x[k][0] - coord[0];
                delr2[1] = x[k][1] - coord[1];
                delr2[2] = x[k][2] - coord[2];
                rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
                if (rsq2 > params[ikparam].cutsq) continue;
                pair->threebody(&params[ijparam],&params[ikparam],&params[ijkparam],rsq1,rsq2,delr1,delr2,fj,fk,eflag,tmp3body);
                threebodyeng+=tmp3body;
            }
        }
        ////////// TOTAL ENERGY/////
        total_energy=threebodyeng+twobodyeng;

        /*
         if (comm->me == 0) {
         printf("in energy(): total_energy= %e, twobodyeng= %e, threebodyeng= %e\n",total_energy,twobodyeng,threebodyeng);
         //printf("end of Stw_GCMC()\n");    // added by Jibao
         }
         */
    } else if (strcmp(pair_style,"hybrid/overlay") != 0) {
        for (int j = 0; j < nall; j++) { // from original lammps; commented by Jibao

            if (i == j) continue;
            if (mode == MOLECULE)
                if (imolecule == molecule[j]) continue;

            delx = coord[0] - x[j][0];
            dely = coord[1] - x[j][1];
            delz = coord[2] - x[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            jtype = type[j];

            if (rsq < cutsq[itype][jtype])
                total_energy +=
                pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
        }   // from original lammps; commented by Jibao
    } else {
        error->all(FLERR,"fix gckmc does not currently support hybrid/overlay pair style");
    }








    /*
     else if (strcmp(pair_style,"hybrid/overlay") != 0) {
     if (!strcmp(pair_style,"sw") || !strcmp(pair_style,"sw0") || !strcmp(pair_style,"sw/omp") || !strcmp(pair_style,"sw0/omp")) {

     } else {

     }
     } else {

     }
     */

    return total_energy;
}

/* ----------------------------------------------------------------------
   compute the energy of the given gas molecule in its current position
   sum across all procs that own atoms of the given molecule
------------------------------------------------------------------------- */

double FixGCkMC::molecule_energy(tagint gas_molecule_id)
{
  double mol_energy = 0.0;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->molecule[i] == gas_molecule_id) {
      mol_energy += energy(i,atom->type[i],gas_molecule_id,atom->x[i]);
    }

  double mol_energy_sum = 0.0;
  MPI_Allreduce(&mol_energy,&mol_energy_sum,1,MPI_DOUBLE,MPI_SUM,world);

  return mol_energy_sum;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixGCkMC::energy_full()
{
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build();
  int eflag = 1;
  int vflag = 0;

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  if (modify->n_post_force) modify->post_force(vflag);
  if (modify->n_end_of_step) modify->end_of_step();

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixGCkMC::pick_random_gas_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int> (ngas*random_equal->uniform());
  if ((iwhichglobal >= ngas_before) &&
      (iwhichglobal < ngas_before + ngas_local)) {
    int iwhichlocal = iwhichglobal - ngas_before;
    i = local_gas_list[iwhichlocal];
  }
  //printf("iwhichglobal = %i, ngas_before = %i, ngas_local = %i\n", iwhichglobal, ngas_before, ngas_local);
  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixGCkMC::pick_random_reactive_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int> (nreact*random_equal->uniform());
  if ((iwhichglobal >= nreact_before) &&
      (iwhichglobal < nreact_before + nreact_local)) {
    int iwhichlocal = iwhichglobal - nreact_before;
    i = local_react_list[iwhichlocal];
  }
  //printf("iwhichglobal = %i, ngas_before = %i, ngas_local = %i\n", iwhichglobal, ngas_before, ngas_local);

  return i;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixGCkMC::pick_random_product_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int> (nprod*random_equal->uniform());
  if ((iwhichglobal >= nprod_before) &&
      (iwhichglobal < nprod_before + nprod_local)) {
    int iwhichlocal = iwhichglobal - nprod_before;
    i = local_prod_list[iwhichlocal];
  }
  //printf("iwhichglobal = %i, ngas_before = %i, ngas_local = %i\n", iwhichglobal, ngas_before, ngas_local);
  return i;
}


/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

tagint FixGCkMC::pick_random_gas_molecule()
{
  int iwhichglobal = static_cast<int> (ngas*random_equal->uniform());
  tagint gas_molecule_id = 0;
  if ((iwhichglobal >= ngas_before) &&
      (iwhichglobal < ngas_before + ngas_local)) {
    int iwhichlocal = iwhichglobal - ngas_before;
    int i = local_gas_list[iwhichlocal];
    gas_molecule_id = atom->molecule[i];
  }

  tagint gas_molecule_id_all = 0;
  MPI_Allreduce(&gas_molecule_id,&gas_molecule_id_all,1,
                MPI_LMP_TAGINT,MPI_MAX,world);

  return gas_molecule_id_all;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGCkMC::toggle_intramolecular(int i)
{
  if (atom->avec->bonds_allow)
    for (int m = 0; m < atom->num_bond[i]; m++)
      atom->bond_type[i][m] = -atom->bond_type[i][m];

  if (atom->avec->angles_allow)
    for (int m = 0; m < atom->num_angle[i]; m++)
      atom->angle_type[i][m] = -atom->angle_type[i][m];

  if (atom->avec->dihedrals_allow)
    for (int m = 0; m < atom->num_dihedral[i]; m++)
      atom->dihedral_type[i][m] = -atom->dihedral_type[i][m];

  if (atom->avec->impropers_allow)
    for (int m = 0; m < atom->num_improper[i]; m++)
      atom->improper_type[i][m] = -atom->improper_type[i][m];
}

/* ----------------------------------------------------------------------
   update the list of gas atoms
------------------------------------------------------------------------- */
//Esteban: asegurarse de que actualize correctamente luego de una reaccion

void FixGCkMC::update_gas_atoms_list()
{
//printf("Begin of FixGCkMC::update_gas_atoms_list()\n");
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  if (nlocal > gcmc_nmax) {
    memory->sfree(local_gas_list);
    gcmc_nmax = atom->nmax;
    local_gas_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "GCMC:local_gas_list");
    memory->sfree(local_react_list);
    gcmc_nmax = atom->nmax;
    local_react_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "GCMC:local_react_list");
    memory->sfree(local_prod_list);
    gcmc_nmax = atom->nmax;
    local_prod_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "GCMC:local_prod_list");
  }

  ngas_local = 0;
 //printf("End of FixGCkMC::update_gas_atoms_list()\n");

}

/* ----------------------------------------------------------------------
   update the list of reactive atoms
------------------------------------------------------------------------- */

void FixGCkMC::update_reactive_atoms_list()
{
 //if (comm->me == 0) printf("Begin of FixGCkMC::update_reactive_atoms_list()\n");
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

//    if (nlocal > gcmc_nmax) {
//    printf("Hasta aca\n");
//    memory->sfree(local_react_list);
//    gcmc_nmax = atom->nmax;
//    local_react_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
//     "GCMC:local_react_list");
//  }

  nreact_local = 0;

    int *type = atom->type; // added by Jibao

  if (regionflag) {

    if (mode == MOLECULE) {

      tagint maxmol = 0;
      for (int i = 0; i < nlocal; i++) maxmol = MAX(maxmol,molecule[i]);
      tagint maxmol_all;
      MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
      double comx[maxmol_all];
      double comy[maxmol_all];
      double comz[maxmol_all];
      for (int imolecule = 0; imolecule < maxmol_all; imolecule++) {
        for (int i = 0; i < nlocal; i++) {
          if (molecule[i] == imolecule) {
            mask[i] |= molecule_group_bit;
          } else {
            mask[i] &= molecule_group_inversebit;
          }
        }
        double com[3];
        com[0] = com[1] = com[2] = 0.0;
        group->xcm(molecule_group,gas_mass,com);
        comx[imolecule] = com[0];
        comy[imolecule] = com[1];
        comz[imolecule] = com[2];
      }

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          if (domain->regions[iregion]->match(comx[molecule[i]],
             comy[molecule[i]],comz[molecule[i]]) == 1) {
            local_gas_list[ngas_local] = i;
            ngas_local++;
          }
        }
      }

    } else { //Esteban: modificado para trabajar con el numero de reactivos
      for (int i = 0; i < nlocal; i++) {
          if ((mask[i] & groupbit) && (type[i] == reactive_type)) {  // Modified by Jibao
        //if (mask[i] & groupbit) { // commented out by Jibao
          if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_react_list[nreact_local] = i;
            nreact_local++;
          }
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == reactive_type)) { // Modified by Jibao
        //if (type[i] == reactive_type) {
      //if (mask[i] & groupbit) {   // commented out by Jibao
        local_react_list[nreact_local] = i;
        nreact_local++;
      }
    }
  }

  MPI_Allreduce(&nreact_local,&nreact,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nreact_local,&nreact_before,1,MPI_INT,MPI_SUM,world);
  nreact_before -= nreact_local;
  //printf("proc=%i, nlocal=%i, nreact_local=%i, nreact_before=%i\n", comm->me, nlocal, nreact_local, nreact_before);
}

/* ----------------------------------------------------------------------
   update the list of product atoms
------------------------------------------------------------------------- */

void FixGCkMC::update_product_atoms_list()
{
 //if (comm->me == 0) printf("Begin of FixGCkMC::update_product_atoms_list()\n");
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

//  if (nlocal > gcmc_nmax) {
//    memory->sfree(local_prod_list);
//    gcmc_nmax = atom->nmax;
//    local_prod_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
//     "GCMC:local_prod_list");
//  }

  nprod_local = 0;

    int *type = atom->type; // added by Jibao

  if (regionflag) {

    if (mode == MOLECULE) {

      tagint maxmol = 0;
      for (int i = 0; i < nlocal; i++) maxmol = MAX(maxmol,molecule[i]);
      tagint maxmol_all;
      MPI_Allreduce(&maxmol,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
      double comx[maxmol_all];
      double comy[maxmol_all];
      double comz[maxmol_all];
      for (int imolecule = 0; imolecule < maxmol_all; imolecule++) {
        for (int i = 0; i < nlocal; i++) {
          if (molecule[i] == imolecule) {
            mask[i] |= molecule_group_bit;
          } else {
            mask[i] &= molecule_group_inversebit;
          }
        }
        double com[3];
        com[0] = com[1] = com[2] = 0.0;
        group->xcm(molecule_group,gas_mass,com);
        comx[imolecule] = com[0];
        comy[imolecule] = com[1];
        comz[imolecule] = com[2];
      }

      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          if (domain->regions[iregion]->match(comx[molecule[i]],
             comy[molecule[i]],comz[molecule[i]]) == 1) {
            local_gas_list[ngas_local] = i;
            ngas_local++;
          }
        }
      }

    } else { //Esteban: modificado para trabajar con el numero de productos
      for (int i = 0; i < nlocal; i++) {
          if ((mask[i] & groupbit) && (type[i] == product_type)) {  // Modified by Jibao
        //if (mask[i] & groupbit) { // commented out by Jibao
          if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_prod_list[nprod_local] = i;
            nprod_local++;
          }
        }
      }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == product_type)) { // Modified by Jibao
        //if (type[i] == reactive_type) {
      //if (mask[i] & groupbit) {   // commented out by Jibao
        local_prod_list[nprod_local] = i;
        nprod_local++;
      }
    }
  }

  MPI_Allreduce(&nprod_local,&nprod,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nprod_local,&nprod_before,1,MPI_INT,MPI_SUM,world);
  nprod_before -= nprod_local;
  //printf("proc=%i, nlocal=%i, nprod_local=%i, nprod_before=%i\n", comm->me, nlocal, nprod_local, nprod_before);
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double FixGCkMC::compute_vector(int n)
{
  if (n == 0) return ntranslation_attempts;
  if (n == 1) return ntranslation_successes;
  if (n == 2) return ninsertion_attempts;
  if (n == 3) return ninsertion_successes;
  if (n == 4) return ndeletion_attempts;
  if (n == 5) return ndeletion_successes;
  if (n == 6) return nrotation_attempts;
  if (n == 7) return nrotation_successes;
    if (n == 8) return energyout;   // added by Jibao
  if (n == 9) return nfreaction_attempts;
  if (n == 10) return nfreaction_successes;
  if (n == 11) return nbreaction_attempts;
  if (n == 12) return nbreaction_successes;

  return 0.0;
}
//Esteban:Agregar los nuevos eventos

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixGCkMC::memory_usage()
{
  double bytes = gcmc_nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixGCkMC::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = random_equal->state();
  list[n++] = random_unequal->state();
  list[n++] = next_reneighbor;
  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixGCkMC::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random_equal->reset(seed);

  seed = static_cast<int> (list[n++]);
  random_unequal->reset(seed);

  next_reneighbor = static_cast<int> (list[n++]);
}
