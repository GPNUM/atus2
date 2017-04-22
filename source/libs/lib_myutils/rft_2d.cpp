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
#include "rft_2d.h"
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace Fourier
{
  /**
   * \brief rft_2d Constructor
   *
   * Real Fourier Transformation in 2 dimension.
   *
   * @param header Header information to construct rft object
   * @param in real input array containing real space data to be transformed
   * @param out complex output array containing transformed data in fourier space
   */
  rft_2d::rft_2d( const generic_header &header, bool b, bool f, Fourier::TYPE t ) : cft_base( header, b, f, t )
  {
    m_forwardPlan  = fftw_plan_dft_r2c_2d( m_dim_x, m_dim_y, m_in_real, m_out, FFTW_ESTIMATE );
    m_backwardPlan = fftw_plan_dft_c2r_2d( m_dim_x, m_dim_y, m_out, m_in_real, FFTW_ESTIMATE );
  }

  /**
   * \brief Performs Fourier Transformation
   *
   * Performs a Fourier Transformation on members in cft_2d object.
   * Forward FT (isign = -1) transforms data in m_in --> m_out
   * Backward FT (isign = 1) transforms data in m_out --> m_in
   *
   * @param isign Whether to perform forward [-1] or backward [1]
   fourier transformation
  */
  void rft_2d::ft( int isign )
  {
    m_isign = isign;
    if ( abs(isign) != 1 ) return;

    if ( isign == -1 )
    {
      fftw_execute( m_forwardPlan );
      scale( m_dx*m_dy/(2*M_PI) );
    }
    else
    {
      fftw_execute( m_backwardPlan );
      scale( m_dkx*m_dky/(2*M_PI) );
    }
  }

  CPoint<2> rft_2d::Get_k(const int64_t l)
  {
    CPoint<2> retval;
    int64_t i = l / m_red_dim;
    int64_t j = l - i*m_red_dim;
    retval[0] = m_dkx*double((i+m_shift_x)%m_dim_x-m_shift_x);
    retval[1] = m_dky*double(j%(m_red_dim-1));
    return retval;
  }

  CPoint<2> rft_2d::Get_x(const int64_t l)
  {
    CPoint<2> retval;
    int64_t i = l / m_dim_y;
    int64_t j = l - i*m_dim_y;
    retval[0] = m_dx*double(i-m_shift_x);
    retval[1] = m_dy*double(j-m_shift_y);
    return retval;
  }

  /**
   * \brief Calculate 1st derivative with respect to x of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void rft_2d::Diff_x()
  {
    ft(-1);

    CPoint<2> k;
    #pragma omp parallel for private(k)
    for ( int64_t i=0; i<m_dim_fs; i++ )
    {
      k = Get_k(i);

      double tmp = m_out[i][0];
      m_out[i][0] = -m_out[i][1]*k[0];
      m_out[i][1] = tmp*k[0];
    }
    ft(1);
  }

  /**
   * \brief Calculate 2nd derivative with respect to x of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void rft_2d::Diff_xx()
  {
    ft(-1);

    CPoint<2> k;
    #pragma omp parallel for private(k)
    for ( int64_t i=0; i<m_dim_fs; i++ )
    {
      k = Get_k(i);

      double tmp = -k[0]*k[0];
      m_out[i][0] *= tmp;
      m_out[i][1] *= tmp;
    }
    ft(1);
  }

  /**
   * \brief Calculate 1st derivative with respect to y of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void rft_2d::Diff_y()
  {
    ft(-1);

    CPoint<2> k;
    #pragma omp parallel for private(k)
    for ( int64_t i=0; i<m_dim_fs; i++ )
    {
      k = Get_k(i);

      double tmp = m_out[i][0];
      m_out[i][0] = -m_out[i][1]*k[1];
      m_out[i][1] = tmp*k[1];
    }
    ft(1);
  }

  /**
   * \brief Calculate 2nd derivative with respect to y of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void rft_2d::Diff_yy()
  {
    ft(-1);

    CPoint<2> k;
    #pragma omp parallel for private(k)
    for ( int64_t i=0; i<m_dim_fs; i++ )
    {
      k = Get_k(i);

      double tmp = -k[1]*k[1];
      m_out[i][0] *= tmp;
      m_out[i][1] *= tmp;
    }
    ft(1);
  }

  /**
   * \brief Applies Laplace Operator on data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void rft_2d::Laplace()
  {
    ft(-1);

    CPoint<2> k;
    #pragma omp parallel for private(k)
    for ( int64_t i=0; i<m_dim_fs; i++ )
    {
      k = Get_k(i);

      double tmp = -(k*k);
      m_out[i][0] *= tmp;
      m_out[i][1] *= tmp;
    }
    ft(1);
  }

  /**
   * \brief Scale data by factor
   *
   * @param faktor Scaling factor
   */
  void rft_2d::scale( const double faktor )
  {
    if ( m_isign == 1 )
    {
      #pragma omp parallel for
      for ( int64_t l=0; l<m_dim; l++ )
      {
        m_in_real[l] *= faktor;
      }
    }
    else
    {
      #pragma omp parallel for
      for ( int64_t l=0; l<m_dim_fs; l++ )
      {
        m_out[l][0] *= faktor;
        m_out[l][1] *= faktor;
      }
    }
  }
}
