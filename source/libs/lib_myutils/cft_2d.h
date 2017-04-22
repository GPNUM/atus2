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


#ifndef CFT_2D_H
#define CFT_2D_H

#include "cft_base.h"

namespace Fourier
{
  /// Class for Fourier transform in two dimensions with complex valued data
  class cft_2d : public Fourier::cft_base<2>
  {
  public:
    cft_2d( const generic_header&, bool=true, bool=false );

    void ft( int isign ); // -1 (forward) oder +1 (backward)

    int Get_dimX() {return m_dim_x;}
    int Get_dimY() {return m_dim_y;}

    void Diff_x();
    void Diff_y();

    void Diff_xx();
    void Diff_yy();

    void Laplace();

    CPoint<2> Get_k(const int64_t) final;
    CPoint<2> Get_x(const int64_t) final;
  protected:
    void fix( fftw_complex* data, double sx, double sy );
    void scale( fftw_complex* data, double sx, double sy );

    void Get_k( const int i,const int j, double & k_x, double & k_y );
    void Get_k( const int i,const int j, int & t_i, int & t_j, double & k_x, double & k_y );
    void Get_kx( const int i, int & t_i, double & k_x );
    void Get_ky( const int j, int & t_j, double & k_y );
    double Get_X(const int i) {return m_dx*(i-m_shift_x);}
    double Get_Y(const int j) {return m_dy*(j-m_shift_y);}
    double Get_kx( const int i );
    double Get_ky( const int j );
  };
}
#endif
