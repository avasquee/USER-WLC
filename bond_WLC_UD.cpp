/* ------------------------------------------------------------------------- */

#include "bond_WLC_UD.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondWLCUD::BondWLCUD(LAMMPS *lmp) : Bond(lmp)
{
  reinitflag = 1;
}

/* ---------------------------------------------------------------------- */

BondWLCUD::~BondWLCUD()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
  }
}

/* ---------------------------------------------------------------------- */

void BondWLCUD::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r, b1,laux1,laux2,laux3;

  ebond = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = (delx*delx + dely*dely + delz*delz)/q0[type]/q0[type];
    r = sqrt(rsq);
    b1 = 1 - rsq;
    laux1=aux1[type];
    laux2=aux2[type];
    laux3=aux3[type];
   
    // force & energy

    if (r > 0.0) fbond = - ( 1.0/b1/b1 - laux1/b1 + laux2 + laux3*b1 )/q0[type]/r0[type];
    else fbond = 0.0;

    if (eflag) ebond = (-1/b1  - laux1*log(b1) - laux2*rsq + laux3*rsq*(1+b1)/2.0)/2.0/r0[type];
    
    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondWLCUD::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");
  memory->create(q0,n+1,"bond:q0");
  memory->create(aux1,n+1,"bond:aux1");
  memory->create(aux2,n+1,"bond:aux2");
  memory->create(aux3,n+1,"bond:aux3");
  
  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondWLCUD::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double r0_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    q0[i] =  k_one*r0_one;
    aux1[i] = 7.0/k_one;
    aux2[i] = 3.0/32.0 - 3.0/4.0/k_one - 6.0/k_one/k_one;
    aux3[i] = ( 13.0/32.0 + 0.81720/k_one - 14.790/k_one/k_one )/( 1.0 - 4.2250/k_one + 4.870/k_one/k_one );
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondWLCUD::equilibrium_distance(int i)
{
  return q0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondWLCUD::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondWLCUD::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k[1],sizeof(double),atom->nbondtypes,fp,NULL,error);
    utils::sfread(FLERR,&r0[1],sizeof(double),atom->nbondtypes,fp,NULL,error);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondWLCUD::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],r0[i]);
}

/* ---------------------------------------------------------------------- */

double BondWLCUD::single(int type, double rsq, int /*i*/, int /*j*/,
                        double &fforce)
{
  double r = sqrt(rsq);
  double b1 = 1 - rsq;
  double laux1=aux1[type];
  double laux2=aux2[type];
  double  laux3=aux3[type];

  fforce = 0;
  if (r > 0.0) fforce = ( 1.0/b1/b1 - laux1/b1 + laux2 + laux3*b1 )/q0[type]/r0[type];
  return (-2/b1 + 2*laux2*b1 + laux3*b1*b1 - 2*laux1*log(r - 1))/q0[type]/r0[type]/4;
}

/* ----------------------------------------------------------------------
    Return ptr to internal members upon request.
------------------------------------------------------------------------ */
void *BondWLCUD::extract( char *str, int &dim )
{
  dim = 1;
  if (strcmp(str,"kappa")==0) return (void*) k;
  if (strcmp(str,"r0")==0) return (void*) r0;
  return NULL;
}


