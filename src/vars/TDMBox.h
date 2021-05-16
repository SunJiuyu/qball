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
//
////////////////////////////////////////////////////////////////////////////////
//
// TDMBox.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef TDMBOX_H
#define TDMBOX_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include <qball/Sample.h>

class TDMBox : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "TDMBox"; };
  
  int set ( int argc, char **argv ) {
    if ( argc == 2 ) {
      int v = atof(argv[1]);
      s->ctrl.TDM_nr_plot = v;

      cout << " <WARNING> TDMBox is set only within unit cell </WARNING>" << endl;
    }

    else if ( argc == 3 ) { 
      int v = atof(argv[1]);
      s->ctrl.TDM_nr_plot = v;
      if(v<1) {
        ui->error("TDM Box [nr] must be larger than 1");
        return 1;
      }


      v = atof(argv[2]);
      s->ctrl.TDM_nR_plot = v;
      if(v<1) {
        ui->error("TDM Box [nR] must be larger than 1");
        return 1;
      }
    }


    else {
      ui->error("TDM Box must be [nr] points in unit cell and [nR] cells");
      return 1;
    }

    return 0;
  }
  
  string print (void) const
  {
    ostringstream st;
    st.setf(ios::left,ios::adjustfield);
    st << setw(10) << name() << " = ";
    st.setf(ios::right,ios::adjustfield);
    st << setw(10) << s->ctrl.TDM_nr_plot;
    st << " in unit cell, and " << setw(10) <<  s->ctrl.TDM_nR_plot << " cells";
    return st.str();
  }
  
  TDMBox(Sample *sample) : s(sample) { s->ctrl.TDM_nr_plot = 2; s->ctrl.TDM_nR_plot = 0 ;}

};

#endif

// Local Variables:
// mode: c++
// End:
