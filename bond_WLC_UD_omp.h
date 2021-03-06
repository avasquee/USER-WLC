/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(WLC/UD/omp,BondWLCUDOMP)

#else

#ifndef LMP_BOND_WLC_UD_OMP_H
#define LMP_BOND_WLC_UD_OMP_H

#include "bond_WLC_UD.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class BondWLCUDOMP : public BondWLCUD, public ThrOMP {

 public:
  BondWLCUDOMP(class LAMMPS *lmp);
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval(int ifrom, int ito, ThrData * const thr);
};

}

#endif
#endif
