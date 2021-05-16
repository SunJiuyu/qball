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
// TDTDM.cc
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TDNEX.cc, time-dependent number of excited electrons -- Jiuyu Jan/2021 $

#include <vector>
#include <valarray>
#include "TDTDM.h"
#include <math/d3vector.h>
#include "Basis.h"
#include <qball/UnitCell.h>
#include "SlaterDet.h"
#include <fstream>
#include "isodate.h"
#include "release.h"
#include "Species.h"
#include <iomanip>
#include "Base64Transcoder.h"
using namespace std;
//#include "NonLocalPotential.h"

TDNEX::TDNEX(const Sample& s, const Wavefunction & wf, VectorPotential * vparg): s_(s), wf_(wf), wf00_(wf), vp(vparg){
  	// initial ground state TDM at t0 in length gauge

  	notGS = false;

   // check initial occ
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for (int ikp=0; ikp<wf_.nkptloc(); ikp++) {

        //double wt = wf_.weight(wf_.kptloc(ikp))
        
        update_projM( *wf_.sd(ispin, ikp) );
      }
    }
  }
  //wf_.wfcontext()->dsum('r',3,1,&psp_current_[0],3);
  	notGS = true;
}


void TDTDM::get_TDNEX(Wavefunction & dwf, bool update_TDM){

	if(!update_TDM) return;

  for ( int ispin = 0; ispin < dwf.nspin(); ispin++ ) {
    if (dwf.spinactive(ispin)) {
      for (int ikp=0; ikp<dwf.nkptloc(); ikp++) {

        double wt = dwf.weight(dwf.kptloc(ikp));
        
        update_projM( *dwf.sd(ispin, ikp) ,wt);
      }
    }
  }

  TranDenMat.resize(nr_plot*nR_plot);
    //cout << "TDTDM good initial" << endl;
	//dwf.wfcontext()->dsum('r',3,1,&DenMat[0],3);
  for (int irR1=0; irR1 < nr_plot*nR_plot; irR1++){
  	TranDenMat[irR1].resize(nr_plot*nR_plot);
  	if ( wf_.context().onpe0() )   cout << "TDTDM: "  << endl;
    for(int irR2=0; irR2 < nr_plot*nR_plot; irR2++){
  		TranDenMat[irR1][irR2] = norm(DenMat[irR1][irR2] - DenMat0[irR1][irR2]);
  		  if ( wf_.context().onpe0() ) {  cout << "TDTDM: " << position_[irR1] <<  " " << position_[irR2] << " " << TranDenMat[irR1][irR2] << endl;}
  	}
  }

}


void TDTDM::update_projM( SlaterDet& sd ){

  complex<double> cZERO = complex<double>(0.0,0.0);
  complex<double> cONE = complex<double>(1.0,0.0);  

          const int nst = sd(ispin,ikp)->nst();
          const int nb = (sd(ispin,ikp)->c()).nb();
          const int nprow = (sd(ispin,ikp)->context()).nprow();
          int mb = nst/nprow + (nst%nprow > 0 ? 1 : 0);


        ComplexMatrix porj_mat(wf_.sd(ispin,ikp)->context(),(wf_.sd(ispin,ikp)->c()).n(),(wf.sd(ispin,ikp)->c()).n(),mb,nb);


    porj_mat.gemm('c','n',1.0,(wf00_.sd(ispin,ikp))->c(),sd(ispin,ikp)->c(),0.0);


  	const ComplexMatrix& c = sd.c();
  	const int nstloc = sd.nstloc();
  	const complex<double>* p = c.cvalptr();
  	const int mloc = c.mloc();
  	const int nbands = c.nloc();
  	const int ngwloc = basis_.localsize();
  	const double * const occ = sd.occ_ptr();



// inner loop for m
    for ( int lj_m=0; lj_m < c.nblocks(); lj_m++ )
      {
        for ( int jj_m=0; jj_m < c.nbs(lj_m); jj_m++ )
        {
          // global state index
          const int mglobal = c.j(lj_m,jj_m);
          const int morig = lj_m*c.nb()+jj_m;


          for ( int ig = 0; ig < ngwloc; ig++ )
          {
            const complex<double> pval = p[ig+norig*mloc];
            const complex<double> pval_m = p_m[ig+morig*mloc_m];             
    		    beta_matrix = conj(pval)*pval_m;                  
          }
	//	    MPI_Barrier(MPI_COMM_WORLD);
        }
      }
//    MPI_Barrier(MPI_COMM_WORLD);
      if (nstloc>0){
        MPI_Comm basis_comm = basis_.context().comm();
//      MPI_Barrier(basis_comm);
        int psi_size = nbands;
        MPI_Allreduce(&dpsi_ik[0], &dpsi_global[0], psi_size,   MPI_DOUBLE_COMPLEX,  MPI_SUM,basis_comm);
        psi_ik_r[inr] = dpsi_global;
      }


}

