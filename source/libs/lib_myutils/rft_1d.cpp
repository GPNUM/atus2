//
// ATUS2 - The ATUS2 package is atom interferometer Toolbox developed at ZARM
// (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
// founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
// 50WM0942, 50WM1042, 50WM1342.
// Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
//
// This file is part of ATUS2.
//
// ATUS2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ATUS2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ATUS2.  If not, see <http://www.gnu.org/licenses/>.
//

#include "rft_1d.h"

namespace Fourier
{
  /**
   * \brief rft_1d Constructor
   *
   * Real Fourier Transformation in 1 dimension.
   *
   * @param header Header information to construct rft object
   */
  rft_1d::rft_1d( const generic_header &header, bool b, bool f, Fourier::TYPE t ) : cft_base( header, b, f, t )
  {
    m_forwardPlan  = fftw_plan_dft_r2c_1d( m_dim, m_in_real, m_out, FFTW_ESTIMATE );
    m_backwardPlan = fftw_plan_dft_c2r_1d( m_dim, m_out, m_in_real, FFTW_ESTIMATE );
  }

  /**
   * \brief Perform Fourier Transformation
   *
   * Performs a Fourier Transformation on members in cft_2d object.
   * Forward FT (isign = -1) transforms data in m_in_real --> m_out
   * Backward FT (isign = 1) transforms data in m_out --> m_in_real
   *
   * @param isign Whether to perform forward [-1] or backward [1]
   fourier transformation
  */
  void rft_1d::ft( int isign )
  {
    int i;
    double faktor;

    m_isign = isign;
    if ( abs(isign) != 1 ) return;
    if ( isign == -1 )
    {
      faktor = m_dx / sqrt(2.0*M_PI);
      fftw_execute( m_forwardPlan );

      for ( i=1; i<m_red_dim; i +=2 )
      {
        m_out[i][0] *= -1;
        m_out[i][1] *= -1;
      }
      for ( i=0; i<m_red_dim; i++   )
      {
        m_out[i][0] *= faktor;
        m_out[i][1] *= faktor;
      }
    }
    else
    {
      faktor = m_dkx / sqrt(2.0*M_PI);

      for ( i=0; i<m_red_dim; i++   )
      {
        m_out[i][0] *= faktor;
        m_out[i][1] *= faktor;
      };
      for ( i=1; i<m_red_dim; i +=2 )
      {
        m_out[i][0] *= -1;
        m_out[i][1] *= -1;
      };

      fftw_execute( m_backwardPlan );
    }
  }

  CPoint<1> rft_1d::Get_k(const int64_t l)
  {
    CPoint<1> retval;
    retval[0] = m_dkx*double(l%(m_red_dim-1));
    return retval;
  }

  CPoint<1> rft_1d::Get_x(const int64_t l)
  {
    CPoint<1> retval;
    retval[0] = double(l-m_shift_x)*m_dx;
    return retval;
  }

  /**
   * \brief Calculates 1st derivative in respect to x on data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void rft_1d::Diff_x()
  {
    ft(-1);

    CPoint<1> k;
    #pragma omp parallel for private(k)
    for ( int i=0; i<m_dim_fs; i++ )
    {
      k = Get_k( i );
      double tmp = m_out[i][0];
      m_out[i][0] = -m_out[i][1]*k[0];
      m_out[i][1] = tmp*k[0];
    }
    ft(1);
  }

  /**
   * \brief Calculates 2nd derivative in respect to x on data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void rft_1d::Diff_xx()
  {
    ft(-1);

    CPoint<1> k;
    #pragma omp parallel for private(k)
    for ( int i=0; i<m_dim_fs; i++ )
    {
      k = Get_k( i );
      double tmp = -(k*k);
      m_out[i][0] *= tmp;
      m_out[i][1] *= tmp;
    }
    ft(1);
  }
}
