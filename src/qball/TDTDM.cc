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
// $Id: TDTDM.cc, time-dependent transition density matrix-- Jiuyu Aug/2020 $

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

TDTDM::TDTDM(const Sample& s, const Wavefunction & wf, VectorPotential * vparg): s_(s), wf_(wf), vp(vparg), nr_plot(s.ctrl.TDM_nr_plot) ,nR_plot(s.ctrl.TDM_nR_plot){
  	// initial ground state TDM at t0 in length gauge


  	r_plot.resize(nr_plot);
  	double a1 = s_.atoms.cell().a(0).x;
  	double a2 = s_.atoms.cell().a(1).y;
  	double a3 = s_.atoms.cell().a(2).z;
  	for(int inr=0; inr<nr_plot; inr++){
  		r_plot[inr].x = 0.51*a1;
  		r_plot[inr].y = 0.51*a2;
  		r_plot[inr].z = inr/(nr_plot*1.0)*a3;
  	}
  	R_plot.resize(nR_plot);
  	for(int inR=0; inR<nR_plot; inR++){
  		R_plot[inR].x = 0.0; //(inR - nR_plot/2.0)*a1;
  		R_plot[inR].y = 0.0; //(inR - nR_plot/2.0)*a2;
  		R_plot[inR].z = (inR*1.0 - nR_plot/2.0)*a3;
  	}

  	DenMat0.resize(nr_plot*nR_plot);
  	DenMat.resize(nr_plot*nR_plot);
  	position_.resize(nr_plot*nR_plot);

  	for (int inr=0; inr< nr_plot; inr++)
  	{
  		for (int inR=0; inR < nR_plot; inR++)
  		{
  			int irR = inr * nR_plot + inR;
  			DenMat0[irR].resize(nr_plot*nR_plot);
  			DenMat[irR].resize(nr_plot*nR_plot);
  			position_[irR] = r_plot[inr].z + R_plot[inR].z ;
  		}
  	}

  	notGS = false;
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ ) {
    if (wf_.spinactive(ispin)) {
      for (int ikp=0; ikp<wf_.nkptloc(); ikp++) {

        double wt = wf_.weight(wf_.kptloc(ikp));
        
        update_DM( *wf_.sd(ispin, ikp) , wt);
      }
    }
  }
  //wf_.wfcontext()->dsum('r',3,1,&psp_current_[0],3);
  	notGS = true;
}


void TDTDM::get_TDTDM(Wavefunction & dwf, bool update_TDM){

	if(!update_TDM) return;

  for ( int ispin = 0; ispin < dwf.nspin(); ispin++ ) {
    if (dwf.spinactive(ispin)) {
      for (int ikp=0; ikp<dwf.nkptloc(); ikp++) {

        double wt = dwf.weight(dwf.kptloc(ikp));
        
        update_DM( *dwf.sd(ispin, ikp) ,wt);
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


void TDTDM::update_DM( SlaterDet& sd , double wt_ik){

  complex<double> cZERO = complex<double>(0.0,0.0);
  complex<double> cONE = complex<double>(1.0,0.0);  


  const Basis& basis_ = sd.basis();
  const double omega = basis_.cell().volume();

    valarray<valarray<complex<double>>> psi_ik_r;
    valarray<valarray<complex<double>>> psi_ik_rR;
    valarray<complex<double>> dpsi_ik, dpsi_global;
//  	vector<complex<double>> psi_r;
//  	psi_r.resize(nr_plot,cZERO);

// DM(r,r',t) = \sum_G(c*exp^{-ikr}) * \sum_G'(c*exp^{ikr'})
  		//			( \psi(r) ) 	*    conjg( \psi(r') )

//  		const double *kpg   = wf_.sd(ispin, ikp)->basis().kpg_ptr();
//  		const double *kpg2  = wf_.sd(ispin, ikp)->basis().kpg2_ptr();
//  		const double *kpgi  = wf_.sd(ispin, ikp)->basis().kpgi_ptr();
  	const double *kpg_x = basis_.kpgx_ptr(0);
  	const double *kpg_y = basis_.kpgx_ptr(1);
  	const double *kpg_z = basis_.kpgx_ptr(2);
//  		const double *gcc_x = wf_.sd(ispin, ikp)->basis().gx_ptr(0);
//  		const double *gcc_y = wf_.sd(ispin, ikp)->basis().gx_ptr(1);
//  		const double *gcc_z = wf_.sd(ispin, ikp)->basis().gx_ptr(2);
  	if(notGS && vp){
    	kpg_x = vp->get_kpgpax(basis_, 0);
    	kpg_y = vp->get_kpgpax(basis_, 1);
    	kpg_z = vp->get_kpgpax(basis_, 2);
  	}



  	const ComplexMatrix& c = sd.c();
  	const int nstloc = sd.nstloc();
  	const complex<double>* p = c.cvalptr();
  	const int mloc = c.mloc();
  	const int nbands = c.nloc();
  	const int ngwloc = basis_.localsize();
  	const double * const occ = sd.occ_ptr();


  	psi_ik_r.resize(nr_plot);
  	psi_ik_rR.resize(nr_plot*nR_plot);
    double kpG_dot_r;
    complex<double> e_ikpGr;
//    e_ikpGr.resize(inr);


	for(int inr=0;inr<nr_plot; inr++){

	  dpsi_ik.resize(nbands,cZERO);
	  dpsi_global.resize(nbands,cZERO);
	  psi_ik_r[inr].resize(nbands,cZERO);

	  for ( int lj=0; lj < c.nblocks(); lj++ )
      {
        for ( int jj=0; jj < c.nbs(lj); jj++ )
        {
          // global state index
          const int nglobal = c.j(lj,jj);
          const int norig = lj*c.nb()+jj;

          for ( int ig = 0; ig < ngwloc; ig++ )
          {
            const complex<double> pval = p[ig+norig*mloc];             
    		kpG_dot_r =  kpg_x[ig]*r_plot[inr].x+kpg_y[ig]*r_plot[inr].y+kpg_z[ig]*r_plot[inr].z;
    		e_ikpGr.real(cos(kpG_dot_r));
    		e_ikpGr.imag(sin(kpG_dot_r));
    		dpsi_ik[nglobal] += pval*e_ikpGr;                  
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

      double k_dot_R;
      complex<double> e_ikR;
      //e_ikR.resize(nR_plot);
      for(int inR=0; inR<nR_plot; inR++){
        k_dot_R =  kpg_x[0]*R_plot[inR].x+kpg_y[0]*R_plot[inR].y+kpg_z[0]*  R_plot[inR].z;        	
        e_ikR.real(cos(k_dot_R));
        e_ikR.imag(sin(k_dot_R));

        int irR = inr * nR_plot + inR;
        psi_ik_rR[irR].resize(nbands,cZERO);
        for(int ib=0; ib< nbands ; ib++) {
    	  psi_ik_rR[irR][ib] = psi_ik_r[inr][ib] * e_ikR;
    	}
      }      

    } //inr


    for (int irR1=0; irR1 < nr_plot*nR_plot; irR1++){
      for(int irR2=0; irR2 < nr_plot*nR_plot; irR2++){
        for(int ib=0; ib< nbands ; ib++) {
        	if(notGS){
        DenMat[irR1][irR2] += 0.5*occ[ib]*wt_ik/omega * psi_ik_rR[irR1][ib]* conj(psi_ik_rR[irR2][ib]) ;
    		}
    		else{
        DenMat0[irR1][irR2] += 0.5*occ[ib]*wt_ik/omega * psi_ik_rR[irR1][ib]* conj(psi_ik_rR[irR2][ib]) ;    			
    		}
        }
      }
    }

}

