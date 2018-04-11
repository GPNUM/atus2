/* * ATUS2 - The ATUS2 package is atom interferometer Toolbox developed at ZARM
 * (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
 * founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
 * 50WM0942, 50WM1042, 50WM1342.
 * Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
 *
 * This file is part of ATUS2.
 *
 * ATUS2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ATUS2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ATUS2.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "my_structs.h"
#include <vector>
#include <cstdlib>

#ifndef ZERNIKE_H
#define ZERNIKE_H

using namespace std;

class zernike 
{
public:
  zernike( const generic_header& );
  ~zernike();

  void rescale_all_zernike( double scale);
  //void generate_zernike();
  void calc_zernike( const generic_header&, double* field );
  //void print_zernike();
  void set_max(double val);
  void add_zernike( int n, int m, double scaling );
protected:
  void resize_zernike( int x, int y );
  vector< vector<double>> m_zern;
  double *m_field;
  generic_header header;
  int m_no_pts;
};
#endif
