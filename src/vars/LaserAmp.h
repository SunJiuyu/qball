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
// LaserAmp.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef LASERAMP_H 
#define LASERAMP_H 

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include <qball/Sample.h>

class LaserAmp : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "laser_amp"; };
  
  int set ( int argc, char **argv ) {
    if ( argc != 4 ) {
      ui->error("laser_amp must be specified as a vector with 3 values");
      return 1;
    }
    
    D3vector v(atof(argv[1]), atof(argv[2]), atof(argv[3]));
    
    s->ctrl.laser_amp= v;
    return 0;
  }
  
  string print (void) const
  {
    ostringstream st;
    st.setf(ios::left,ios::adjustfield);
    st << setw(10) << name() << " = ";
    st.setf(ios::right,ios::adjustfield);
    st << setw(10) << s->ctrl.laser_amp;
    return st.str();
  }
  
  LaserAmp(Sample *sample) : s(sample) { s->ctrl.laser_amp= D3vector(0.0, 0.0, 0.0); }

};

#endif

// Local Variables:
// mode: c++
// End:
