/* ------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(WLC/UD,BondWLCUD)

#else

#ifndef LMP_BOND_WLC_UD_H
#define LMP_BOND_WLC_UD_H

#include "bond.h"

namespace LAMMPS_NS {

class BondWLCUD : public Bond {
 public:
  BondWLCUD(class LAMMPS *);
  virtual ~BondWLCUD();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  virtual void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);
  virtual void *extract(char *, int &);

 protected:
  double *k,*r0,*q0,*aux1,*aux2,*aux3;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
