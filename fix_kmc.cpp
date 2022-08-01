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

#include "fix_kmc.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "compute.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "memory.h"
#include "error.h"
#include "output.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_park.h"
#include "math_const.h"
#include <iostream>
#include <cmath>
#include <cstring>

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{ATOM,MOLECULE};

/* ---------------------------------------------------------------------- */

Fixkmc::Fixkmc(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  list(nullptr)
{
  if (narg < 11) 
    error->all(FLERR,"Incorrect number of fix kmc arguments {}", narg);

  if (atom->molecular == 2)
    error->all(FLERR,"Fix kmc does not (yet) work with atom_style template");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args
  reservoir_temperature = utils::numeric(FLERR,arg[3],false,lmp);
  reactive_type = utils::inumeric(FLERR,arg[4],false,lmp);
  product_type = utils::inumeric(FLERR,arg[5],false,lmp);
  surf_type = utils::inumeric(FLERR,arg[6],false,lmp);
  preexp = utils::numeric(FLERR,arg[7],false,lmp);
  potential = utils::numeric(FLERR,arg[8],false,lmp);
  seed = utils::inumeric(FLERR,arg[9],false,lmp);
  nevery = utils::inumeric(FLERR,arg[10],false,lmp);

  if (seed <= 0) 
    error->all(FLERR,"Illegal fix kmc seed {}", seed);
  if (reservoir_temperature < 0.0) 
    error->all(FLERR,"Illegal fix kmc reservoir temperature {}", reservoir_temperature);

  regionflag = 0;

  // read options from end of input line
  options(narg-11,&arg[11]);

  // random number generator, same for all procs
  random_equal = new RanPark(lmp,seed);

  // random number generator, not the same for all procs
  random_unequal = new RanPark(lmp,seed);

  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi =
    region_zlo = region_zhi = 0.0;
  if (regionflag) {
    if (iregion->bboxflag == 0)
      error->all(FLERR,"Fix kmc region does not support a bounding box");
    if (iregion->dynamic_check())
      error->all(FLERR,"Fix kmc region cannot be dynamic");

    region_xlo = iregion->extent_xlo;
    region_xhi = iregion->extent_xhi;
    region_ylo = iregion->extent_ylo;
    region_yhi = iregion->extent_yhi;
    region_zlo = iregion->extent_zlo;
    region_zhi = iregion->extent_zhi;

    if (region_xlo < domain->boxlo[0] || region_xhi > domain->boxhi[0] ||
        region_ylo < domain->boxlo[1] || region_yhi > domain->boxhi[1] ||
        region_zlo < domain->boxlo[2] || region_zhi > domain->boxhi[2])
      error->all(FLERR,"Fix kmc region extends outside simulation box");

  }

  if (mode == ATOM) natoms_per_molecule = 1;
  else error->all(FLERR,"Fix kmc region does not support molecules");

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep+1;

  // zero out counters
  nfreaction_attempts = 0.0;
  nfreaction_successes = 0.0;
  nbreaction_attempts = 0.0;
  nbreaction_successes = 0.0;

  gcmc_nmax = 0;
  local_gas_list = NULL;
  local_react_list = NULL;
  local_prod_list = NULL;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void Fixkmc::options(int narg, char **arg)
{
  if (narg < 0) 
    utils::missing_cmd_args(FLERR, "fix kmc", error);

  // defaults

  mode = ATOM;
  regionflag = 0;
  iregion = nullptr;

  int iarg = 0;
  while (iarg < narg) {
  if (strcmp(arg[iarg],"mol") == 0) {
      error->all(FLERR,"Fix kmc does not work with molecules (yet)!");
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) 
        utils::missing_cmd_args(FLERR, "fix kmc", error);
      iregion = domain->get_region_by_id(arg[iarg+1]);
      if (iregion == nullptr)
        error->all(FLERR,"Region ID for fix kmc does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix kmc command {}", arg[iarg]);
  }
}

/* ---------------------------------------------------------------------- */

Fixkmc::~Fixkmc()
{
  if (regionflag) delete [] idregion;
  delete random_equal;
  delete random_unequal;

  memory->destroy(local_gas_list);
  memory->destroy(local_react_list);
  memory->destroy(local_prod_list);

}

/* ---------------------------------------------------------------------- */

int Fixkmc::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}


/* ---------------------------------------------------------------------- */

void Fixkmc::init()
{
  triclinic = domain->triclinic;

  int *type = atom->type;
  if (mode == ATOM) {
    if (product_type <= 0 || product_type > atom->ntypes)
      error->all(FLERR,"Invalid product atom type in fix kmc command {}", product_type);
    if (reactive_type <= 0 || reactive_type > atom->ntypes)
      error->all(FLERR,"Invalid reactive atom type in fix kmc command {}", reactive_type);
  }

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix kmc in a 2d simulation");

  // create a new group for interaction exclusions

  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] == groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot do kmc on atoms in atom_modify first group");
  }

  beta = 1.0/(1.38065e-23*reservoir_temperature);

  kfreact = (float)nevery*preexp*exp(0.5*1*1.6022e-19*beta*potential);
  kbreact = (float)nevery*preexp*exp(-0.5*1*1.6022e-19*beta*potential);

  imagezero = ((imageint) IMGMAX << IMG2BITS) |
             ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // construct group bitmask for all new atoms
  groupbitall = 1 | groupbit;

  neighbor->add_request(this,NeighConst::REQ_FULL);
}
/* ---------------------------------------------------------------------- */

void Fixkmc::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   Attempt Kinetic Monte Carlo reactions consisting in changing types
   depending to the distance to a catalyst
------------------------------------------------------------------------- */

void Fixkmc::pre_exchange()
{
  if (next_reneighbor != update->ntimestep) return;
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

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();

  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  update_gas_atoms_list();
  update_product_atoms_list();
  update_reactive_atoms_list();

  attempt_atomic_freaction(nreact);

  next_reneighbor = update->ntimestep + nevery;

}

void Fixkmc::attempt_atomic_freaction(int nreact)
{
  nfreaction_attempts += 1.0;
  int success = 0;

  int *ilist,*jlist,*numneigh,**firstneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz, kvel;
  double rsq,rsq1;

  double **x = atom->x;
  int *type = atom->type;
  int tstep = update->dt;
  int inum, jnum, k;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for(int i; i<nreact; i++)
  {
    // is this atom in this processor
    if ((i >= nreact_before) && (i < nreact_before + nreact_local)){  

      int ilocal= i - nreact_before;

      // j is the actual id of the atom
      int j = local_react_list[ilocal];  

      xtmp = x[j][0];
      ytmp = x[j][1];
      ztmp = x[j][2];

      jlist = firstneigh[j];
      jnum = numneigh[j];
      rsq1 = INT_MAX;

      // go find me the closest surf_type
      for (int jj = 0; jj < jnum; jj++) { 
        k = jlist[jj];
        k &= NEIGHMASK;
        if (type[k] == surf_type)
        {
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < rsq1) rsq1 = rsq;
        }
      }

    rsq1 = sqrt(rsq1);
    kvel = exp(-rsq1)*kfreact;

      if (random_unequal->uniform() <
          1-exp(-kvel*tstep)) {
              type[j] = product_type;
              success += 1;
      }
    }

    // Comunicate to the other processors and remake lists if needed
    int success_all = 0;
    MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);
    if (success_all) {
        update_gas_atoms_list();
        update_reactive_atoms_list();
        update_product_atoms_list();
        nfreaction_successes += 1;

        atom->nghost = 0;
        comm->borders();
        comm->exchange();
    }
  }
}

/* ----------------------------------------------------------------------
   update the list of gas atoms. Esteban: Now just initializes the memory
   for the other lists
------------------------------------------------------------------------- */

void Fixkmc::update_gas_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  if (nlocal > gcmc_nmax) {
    memory->sfree(local_gas_list);
    gcmc_nmax = atom->nmax;
    local_gas_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmc:local_gas_list");
    memory->sfree(local_react_list);
    gcmc_nmax = atom->nmax;
    local_react_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmc:local_react_list");
    memory->sfree(local_prod_list);
    gcmc_nmax = atom->nmax;
    local_prod_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmc:local_prod_list");
  }

  ngas_local = 0;

}

/* ----------------------------------------------------------------------
   update the list of reactive atoms
------------------------------------------------------------------------- */

void Fixkmc::update_reactive_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  nreact_local = 0;

  int *type = atom->type;

  if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == reactive_type)) {
          if (iregion->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_react_list[nreact_local] = i;
            nreact_local++;
          }
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == reactive_type)) {
        local_react_list[nreact_local] = i;
        nreact_local++;
      }
    }
  }

  MPI_Allreduce(&nreact_local,&nreact,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nreact_local,&nreact_before,1,MPI_INT,MPI_SUM,world);
  nreact_before -= nreact_local;
}

/* ----------------------------------------------------------------------
   update the list of product atoms
------------------------------------------------------------------------- */

void Fixkmc::update_product_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  nprod_local = 0;

  int *type = atom->type;

  if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == product_type)) {
          if (iregion->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_prod_list[nprod_local] = i;
            nprod_local++;
          }
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == product_type)) {
        local_prod_list[nprod_local] = i;
        nprod_local++;
      }
    }
  }

  MPI_Allreduce(&nprod_local,&nprod,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nprod_local,&nprod_before,1,MPI_INT,MPI_SUM,world);
  nprod_before -= nprod_local;
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double Fixkmc::compute_vector(int n)
{
  if (n == 0) return nfreaction_attempts;
  if (n == 1) return nfreaction_successes;
  if (n == 2) return nbreaction_attempts;
  if (n == 3) return nbreaction_successes;

  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double Fixkmc::memory_usage()
{
  double bytes = gcmc_nmax * sizeof(int);
  return bytes;
}
