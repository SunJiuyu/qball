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
// $Id: TDTDM.h, time-dependent transition density matrix-- Jiuyu Aug/2020 $

#ifndef TDTDM_H
#define TDTDM_H

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

class TDTDM 
{
  private:

  const Wavefunction& wf_;
  const Sample& s_;
  //const Basis& basis_;

  VectorPotential * vp;

  vector<double> position_;


  bool notGS;


  public:

  int nr_plot;
  int nR_plot;
  valarray<D3vector> r_plot;
  valarray<D3vector> R_plot;  
  
  vector<vector<double> > TranDenMat;

  vector<vector<complex<double>> > DenMat0; 
  vector<vector<complex<double>> > DenMat; 


  TDTDM (const Sample & s, const Wavefunction & wf, VectorPotential * vparg);
  ~TDTDM (){
  }

  void get_TDTDM( Wavefunction & dwf, bool update_TDM);
  void update_DM( SlaterDet& sd, double wt_ik);

};

#endif
