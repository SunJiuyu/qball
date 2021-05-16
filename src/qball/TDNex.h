////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////
//
// TDTDM.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDTDM.h, time-dependent number of excited electrons -- Jiuyu Jan/2021 $

#ifndef TDNEX_H
#define TDNEX_H

#include <vector>
#include <valarray>
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include <math/matrix.h>
#include <math/d3vector.h>
#include "PrintMem.h"
#include <math/blas.h>
using namespace std;

class TDNEX 
{
  private:

  const Wavefunction& wf_;
  const Wavefunction wf00_; // for GS wavefuntion storage -- Jiuyu
  const Sample& s_;
  //const Basis& basis_;

  VectorPotential * vp;


  bool notGS;


  public:

  
  vector<vector<double> > TranDenMat;

  vector<vector<complex<double>> > proj_mat; 
  vector<vector<complex<double>> > num_exc_ele; 


  TDNEX (const Sample & s, const Wavefunction & wf, VectorPotential * vparg);
  ~TDNEX (){
  }

  void get_TDNEX( Wavefunction & dwf, bool update_TDM);
  void update_projM( SlaterDet& sd, double wt_ik); // update project matrix

};

#endif
