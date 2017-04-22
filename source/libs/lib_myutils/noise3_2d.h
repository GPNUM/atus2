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

#include <ostream>
#include "fftw3.h"
#include "cft_2d.h"
#include "cft_1d.h"
#include "my_structs.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

using namespace std;

namespace Fourier
{
  class CNoise2_2D : public Fourier::cft_2d
  {
  public:
    CNoise2_2D( const generic_header& );

    void Do_Noise_Mirror( double sigma_0 = 10, double rho = 0.1, int p = 3 );
    void Do_Noise_Metric( double max_noise, double sigma = 1.0, double exp_x = 2.0,
                          double exp_y = 2.0 );

    void Set_Seed( const double a ) { m_seed = a; };
    void Set_Alphas( const double a, const double b ) { m_alpha_x = a; m_alpha_y = b; };

    double Get_Val_re( const int i, const int j );
    double Get_Val_im( const int i, const int j );

    double Get_Val( const int ij );


  protected:
    double m_seed;
    double m_alpha_x;
    double m_alpha_y;
  };
}
