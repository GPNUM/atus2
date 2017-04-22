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

#include "cft_base.h"

namespace Fourier
{
  /// Class for Fourier transform in two dimensions with real valued data
  class rft_2d : public cft_base<2>
  {
  public:
    rft_2d( const generic_header&, bool=true, bool=true, Fourier::TYPE=Fourier::TYPE::REAL );

    void ft( int isign ); // -1 (forward) oder +1 (backward)
    void Diff_x();
    void Diff_y();

    void Diff_xx();
    void Diff_yy();
    void Laplace();

    CPoint<2> Get_k(const int64_t) final;
    CPoint<2> Get_x(const int64_t) final;
  private:
    void scale( const double );
  };
}
