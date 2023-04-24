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
   Contributing author: Esteban Gadea (U of Buenos Aires)
------------------------------------------------------------------------- */

#include "fix_kmcbi.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "compute.h"
#include "force.h"
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

Fixkmcbi::Fixkmcbi(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  list(nullptr)
{
  if (narg < 11)  //original 11
    error->all(FLERR,"Incorrect number of fix kmcbi arguments {}", narg);

  if (atom->molecular == 2)
    error->all(FLERR,"Fix kmcbi does not (yet) work with atom_style template");

  dynamic_group_allow = 1;

  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  time_depend = 1;

  // required args
  reservoir_temperature = utils::numeric(FLERR,arg[3],false,lmp);
  reactivea_type = utils::inumeric(FLERR,arg[4],false,lmp);
  reactiveb_type = utils::inumeric(FLERR,arg[5],false,lmp);
  productc_type = utils::inumeric(FLERR,arg[6],false,lmp);
  productd_type = utils::inumeric(FLERR,arg[7],false,lmp);
  //surf_type = utils::inumeric(FLERR,arg[6],false,lmp); 
  preexp = utils::numeric(FLERR,arg[8],false,lmp);
  seed = utils::inumeric(FLERR,arg[9],false,lmp);
  nevery = utils::inumeric(FLERR,arg[10],false,lmp);

  if (seed <= 0)
    error->all(FLERR,"Illegal fix kmcbi seed {}", seed);
  if (reservoir_temperature < 0.0)
    error->all(FLERR,"Illegal fix kmcbi reservoir temperature {}", reservoir_temperature);

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
      error->all(FLERR,"Fix kmcbi region does not support a bounding box");
    if (iregion->dynamic_check())
      error->all(FLERR,"Fix kmcbi region cannot be dynamic");

    region_xlo = iregion->extent_xlo;
    region_xhi = iregion->extent_xhi;
    region_ylo = iregion->extent_ylo;
    region_yhi = iregion->extent_yhi;
    region_zlo = iregion->extent_zlo;
    region_zhi = iregion->extent_zhi;

    if (region_xlo < domain->boxlo[0] || region_xhi > domain->boxhi[0] ||
        region_ylo < domain->boxlo[1] || region_yhi > domain->boxhi[1] ||
        region_zlo < domain->boxlo[2] || region_zhi > domain->boxhi[2])
      error->all(FLERR,"Fix kmcbi region extends outside simulation box");

  }
  if (densityflag){
    if (iregion->bboxflag == 0)
      error->all(FLERR,"Fix kmcbi density does not support a bounding box");
    if (iregion->dynamic_check())
      error->all(FLERR,"Fix kmcbi density cannot be dynamic");

    region_xlo = jregion->extent_xlo;
    region_xhi = jregion->extent_xhi;
    region_ylo = jregion->extent_ylo;
    region_yhi = jregion->extent_yhi;
    region_zlo = jregion->extent_zlo;
    region_zhi = jregion->extent_zhi;

    if (region_xlo < domain->boxlo[0] || region_xhi > domain->boxhi[0] ||
        region_ylo < domain->boxlo[1] || region_yhi > domain->boxhi[1] ||
        region_zlo < domain->boxlo[2] || region_zhi > domain->boxhi[2])
      error->all(FLERR,"Fix kmcbi region extends outside simulation box");

    densvol = (region_xhi-region_xlo)*
              (region_yhi-region_ylo)*
              (region_zhi-region_zlo);
  }

  if (mode == ATOM) natoms_per_molecule = 1;
  else error->all(FLERR,"Fix kmcbi region does not support molecules");

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep+1;

  // zero out counters
  nfreaction_attempts = 0.0;
  nfreaction_successes = 0.0;
  nbreaction_attempts = 0.0;
  nbreaction_successes = 0.0;

  gcmc_nmax = 0;
  local_gas_list = NULL;
  local_reacta_list = NULL;
  local_reactb_list = NULL;
  local_prodc_list = NULL;
  local_prodd_list = NULL;
  local_gap_list = NULL;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void Fixkmcbi::options(int narg, char **arg)
{
  if (narg < 0)
    utils::missing_cmd_args(FLERR, "fix kmcbi", error);

  // defaults

  mode = ATOM;
  regionflag = 0;
  densityflag = 0;
  iregion = nullptr;
  jregion = nullptr;

  int iarg = 0;
  while (iarg < narg) {
  if (strcmp(arg[iarg],"mol") == 0) {
      error->all(FLERR,"Fix kmcbi does not work with molecules (yet)!");
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, "fix kmcbi", error);
      iregion = domain->get_region_by_id(arg[iarg+1]);
      if (iregion == nullptr)
        error->all(FLERR,"Region ID for fix kmcbi does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"density") == 0) {
      if (iarg+4 > narg)
        utils::missing_cmd_args(FLERR, "fix kmcbi", error);
      jregion = domain->get_region_by_id(arg[iarg+1]);
      mindens = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      maxdens = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (jregion == nullptr)
        error->all(FLERR,"Region ID for fix kmcbi does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      jidregion = new char[n];
      strcpy(jidregion,arg[iarg+1]);
      densityflag = 1;
      iarg += 4;
    } else error->all(FLERR,"Illegal fix kmcbi command {}", arg[iarg]);
  }
}

/* ---------------------------------------------------------------------- */

Fixkmcbi::~Fixkmcbi()
{
  if (regionflag) delete [] idregion;
  delete random_equal;
  delete random_unequal;

  memory->destroy(local_gas_list);  
  memory->destroy(local_reacta_list);
  memory->destroy(local_reactb_list);
  memory->destroy(local_prodc_list);
  memory->destroy(local_prodd_list);
  memory->destroy(local_gap_list);

}

/* ---------------------------------------------------------------------- */

int Fixkmcbi::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}


/* ---------------------------------------------------------------------- */

void Fixkmcbi::init()
{
  triclinic = domain->triclinic;

  int *type = atom->type;
  if (mode == ATOM) {
    
    if (reactivea_type <= 0 || reactivea_type > atom->ntypes)
      error->all(FLERR,"Invalid reactive atom type in fix kmcbi command {}", reactivea_type);
    if (reactiveb_type <= 0 || reactiveb_type > atom->ntypes)
      error->all(FLERR,"Invalid reactive atom type in fix kmcbi command {}", reactiveb_type);
    if (productc_type <= 0 || productc_type > atom->ntypes)
      error->all(FLERR,"Invalid product atom type in fix kmcbi command {}", productc_type);
    if (productd_type <= 0 || productd_type > atom->ntypes)
      error->all(FLERR,"Invalid product atom type in fix kmcbi command {}", productd_type);
  }

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix kmcbi in a 2d simulation");

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
      error->all(FLERR,"Cannot do kmcbi on atoms in atom_modify first group");
  }

  beta = 1.0/(force->boltz*reservoir_temperature);

  kfreact = (float)nevery * preexp;
  kbreact = (float)nevery * preexp;

  imagezero = ((imageint) IMGMAX << IMG2BITS) |
             ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // construct group bitmask for all new atoms
  groupbitall = 1 | groupbit;

  neighbor->add_request(this,NeighConst::REQ_FULL);
}
/* ---------------------------------------------------------------------- */

void Fixkmcbi::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   Attempt Kinetic Monte Carlo reactions consisting in changing types
   depending to the distance to a catalyst
------------------------------------------------------------------------- */

void Fixkmcbi::pre_exchange()
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
  update_productc_atoms_list();
  update_productd_atoms_list();
  update_reactivea_atoms_list();
  update_reactiveb_atoms_list();
  update_gap_atoms_list();

  attempt_atomic_freaction(nreacta);

  next_reneighbor = update->ntimestep + nevery;

}

void Fixkmcbi::attempt_atomic_freaction(int nreacta)
{
  int success;
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

  for(int i=0; i<nreacta; i++)
  {
    nfreaction_attempts += 1.0;
    success = 0;
    k = 0;
    // is this atom in this processor
    if ((i >= nreacta_before) && (i < nreacta_before + nreacta_local)){

      int ilocal= i - nreacta_before;

      // j is the actual id of the atom
      int j = local_reacta_list[ilocal];

      xtmp = x[j][0];
      ytmp = x[j][1];
      ztmp = x[j][2];

      jlist = firstneigh[j];
      jnum = numneigh[j];
      rsq1 = INT_MAX;

      // go find me all the reactB neighbors
      for (int jj = 0; jj < jnum; jj++) {
        k = jlist[jj];
        k &= NEIGHMASK;
        if (type[k] == reactiveb_type)
        {
          delx = xtmp - x[k][0];
          dely = ytmp - x[k][1];
          delz = ztmp - x[k][2];
          rsq = delx*delx + dely*dely + delz*delz;
          //if (rsq < rsq1) rsq1 = rsq;
          kvel = exp(-sqrt(rsq))*kfreact;
          if (random_unequal->uniform() < 1-exp(-kvel*tstep)) {
            type[j] = productc_type;
            //shift_type(k);
            type[k] = productd_type;
            success += 1;
            printf("freaction:  r: %6.4f \n",
              sqrt(rsq));
          }

        }
      }
      if (success > 1) printf("WARNING: multiple reactions accepted -> Nsuccess = %i\n Nreacta = %i\n Nreactb = %i\n Nprodc = %i\n Nprodd = %i\n",
        success, nreacta, nreactb, nprodc, nprodd);

      //kvel = exp(-sqrt(rsq1))*kfreact;

      //if (random_unequal->uniform() < 1-exp(-kvel*tstep)) {
      //  type[j] = product_type;
      //  success += 1;
      //  printf("freaction:  x: %6.4f y: %6.4f z: %6.4f\n",
      //    xtmp, ytmp, ztmp);
      //}
    }

    // Comunicate to the other processors and remake lists if needed
    int success_all = 0;
    MPI_Allreduce(&success,&success_all,1,MPI_INT,MPI_MAX,world);
    //int k_all = 0;
    //MPI_Allreduce(&k,&k_all,1,MPI_INT,MPI_MAX,world);
    if (success_all) {
      //shift_type(k_all);
      if (densityflag){
        if ((float)ngap / densvol > mindens){
          swap_random_gap_atom();
        }
        if ((float)ngap / densvol > maxdens){
          swap_random_gap_atom();
        }
      }
        update_gas_atoms_list();
        update_reactivea_atoms_list();
        update_reactiveb_atoms_list();
        update_productc_atoms_list();
        update_productd_atoms_list();
        update_gap_atoms_list();
        nfreaction_successes += 1;
        printf("Nreacta = %i\n Nreactb = %i\n Nprodc = %i\n Nprodd = %i\n",
         nreacta, nreactb, nprodc, nprodd);

        atom->nghost = 0;
        comm->borders();
        comm->exchange();
    }
  }
}


void Fixkmcbi::shift_type(int k)
{
  if ((k >= nreacta_before) && (k < nreacta_before + nreacta_local)){
    int *type = atom->type;
    type[k] = productd_type;
  }
}

/* ----------------------------------------------------------------------
   update the list of gas atoms. Esteban: Now just initializes the memory
   for the other lists
------------------------------------------------------------------------- */

void Fixkmcbi::update_gas_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  if (nlocal > gcmc_nmax) {
    memory->sfree(local_gas_list);
    gcmc_nmax = atom->nmax;
    local_gas_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmcbi:local_gas_list");
    memory->sfree(local_reacta_list);
    gcmc_nmax = atom->nmax;
    local_reacta_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmcbi:local_reacta_list");
    memory->sfree(local_reactb_list);
    gcmc_nmax = atom->nmax;
    local_reactb_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmcbi:local_reactb_list");
    memory->sfree(local_prodc_list);
    gcmc_nmax = atom->nmax;
    local_prodc_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmcbi:local_prodc_list");
    memory->sfree(local_prodd_list);
    gcmc_nmax = atom->nmax;
    local_prodd_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmcbi:local_prodd_list");

    memory->sfree(local_gap_list);
    gcmc_nmax = atom->nmax;
    local_gap_list = (int *) memory->smalloc(gcmc_nmax*sizeof(int),
     "kmcbi:local_gap_list");

  }

  ngas_local = 0;

}

/* ----------------------------------------------------------------------
   update the list of reactive atoms
------------------------------------------------------------------------- */

void Fixkmcbi::update_reactivea_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  nreacta_local = 0;

  int *type = atom->type;

  if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == reactivea_type)) {
          if (iregion->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_reacta_list[nreacta_local] = i;
            nreacta_local++;
          }
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == reactivea_type)) {
        local_reacta_list[nreacta_local] = i;
        nreacta_local++;
      }
    }
  }

  MPI_Allreduce(&nreacta_local,&nreacta,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nreacta_local,&nreacta_before,1,MPI_INT,MPI_SUM,world);
  nreacta_before -= nreacta_local;
}

void Fixkmcbi::update_reactiveb_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  nreactb_local = 0;

  int *type = atom->type;

  if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == reactiveb_type)) {
          if (iregion->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_reactb_list[nreactb_local] = i;
            nreactb_local++;
          }
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == reactiveb_type)) {
        local_reactb_list[nreactb_local] = i;
        nreactb_local++;
      }
    }
  }

  MPI_Allreduce(&nreactb_local,&nreactb,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nreactb_local,&nreactb_before,1,MPI_INT,MPI_SUM,world);
  nreactb_before -= nreactb_local;
}

/* ----------------------------------------------------------------------
   update the list of product atoms
------------------------------------------------------------------------- */

void Fixkmcbi::update_productc_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  nprodc_local = 0;

  int *type = atom->type;

  if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == productc_type)) {
          if (iregion->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_prodc_list[nprodc_local] = i;
            nprodc_local++;
          }
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == productc_type)) {
        local_prodc_list[nprodc_local] = i;
        nprodc_local++;
      }
    }
  }

  MPI_Allreduce(&nprodc_local,&nprodc,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nprodc_local,&nprodc_before,1,MPI_INT,MPI_SUM,world);
  nprodc_before -= nprodc_local;
}


void Fixkmcbi::update_productd_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  nprodd_local = 0;

  int *type = atom->type;

  if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == productd_type)) {
          if (iregion->match(x[i][0],x[i][1],x[i][2]) == 1) {
            local_prodd_list[nprodd_local] = i;
            nprodd_local++;
          }
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++) {
        if ((mask[i] & groupbit) && (type[i] == productd_type)) {
        local_prodd_list[nprodd_local] = i;
        nprodd_local++;
      }
    }
  }

  MPI_Allreduce(&nprodd_local,&nprodd,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&nprodd_local,&nprodd_before,1,MPI_INT,MPI_SUM,world);
  nprodd_before -= nprodd_local;
}

void Fixkmcbi::update_gap_atoms_list()
{
  if (densityflag){
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  double **x = atom->x;

  ngap_local = 0;

  int *type = atom->type;


  for (int i = 0; i < nlocal; i++) {
    if ((mask[i] & groupbit) && (type[i] == productc_type)) {
      if (jregion->match(x[i][0],x[i][1],x[i][2]) == 1) {
        local_gap_list[ngap_local] = i;
        ngap_local++;
      }
    }
  }


  MPI_Allreduce(&ngap_local,&ngap,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ngap_local,&ngap_before,1,MPI_INT,MPI_SUM,world);
  ngap_before -= ngap_local;
}
}

void Fixkmcbi::swap_random_gap_atom()
{
  int *type = atom->type;
  int i = -1;
  int iwhichglobal = static_cast<int> (ngap*random_equal->uniform());
  if ((iwhichglobal >= ngap_before) &&
      (iwhichglobal < ngap_before + ngap_local)) {
    int iwhichlocal = iwhichglobal - ngap_before;
    i = local_gap_list[iwhichlocal];
    type[i] = reactivea_type;
  }
}

/* ----------------------------------------------------------------------
  return acceptance ratios
------------------------------------------------------------------------- */

double Fixkmcbi::compute_vector(int n)
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

double Fixkmcbi::memory_usage()
{
  double bytes = gcmc_nmax * sizeof(int);
  return bytes;
}
