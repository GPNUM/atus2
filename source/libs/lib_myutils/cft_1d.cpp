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


#include <cmath>
#include <cstdlib>
#include <cstring>
#include "cft_1d.h"

namespace Fourier
{
  /**
   * \brief cft_1d Constructor
   *
   * Complex Fourier Transformation in 1 dimension.
   *
   * @param header Header information to construct cft object
   * @param b Whether inplace transformation is done
   */
  cft_1d::cft_1d( const generic_header &header, bool b, bool f ) : cft_base( header, b, f )
  {
    m_bfix = true;

    m_forwardPlan  = fftw_plan_dft_1d( m_dim, m_in, m_out, FFTW_FORWARD, FFTW_ESTIMATE );
    m_backwardPlan = fftw_plan_dft_1d( m_dim, m_out, m_in, FFTW_BACKWARD, FFTW_ESTIMATE );

    assert( m_forwardPlan != nullptr );
    assert( m_backwardPlan != nullptr );
  }

  /**
   * \brief Performs Fourier Transformation
   *
   * Performs a Fourier Transformation on members in cft_1d object.
   * Forward FT (isign = -1) transforms data in m_in --> m_out
   * Backward FT (isign = 1) transforms data in m_out --> m_in
   *
   * @param isign Whether to perform forward [-1] or backward [1]
   fourier transformation
  */
  void cft_1d::ft( int isign )
  {
    m_isign = isign;
    if ( abs(isign) != 1 ) return;
    if ( isign == -1 )
    {
      fftw_execute( m_forwardPlan );
      if ( m_bfix ) fix( m_out, m_dx );
      else scale( m_out, m_dx );
    }
    else
    {
      fftw_execute( m_backwardPlan );
      if ( m_bfix ) fix( m_in, m_dkx );
      else scale( m_in, m_dkx );
    }
  }

  /**
  * \brief Get x position
  *
  * @param[in] linear array index
  * @returns CPoint<1> containing the spatial position
  */
  CPoint<1> cft_1d::Get_x( const int64_t l )
  {
    CPoint<1> retval;
    retval[0] = double(l-m_shift_x)*m_dx;
    return retval;
  }

  /**
   * \brief Get k value
   *
   * @param[in] linear array index
   * @returns CPoint<1> containing the spatial frequency
   */
  CPoint<1> cft_1d::Get_k( const int64_t l )
  {
    CPoint<1> retval;
    int64_t t_i;
    if ( !m_bfix )
      t_i = (l+m_shift_x)%m_dim;
    else
      t_i = l;
    retval[0] = double(t_i-m_shift_x)*m_dkx;
    return retval;
  }

  /**
  * \brief Get k value
  *
  * @param[in] linear array index
  * @returns CPoint<1> containing the frequency
  */
  void cft_1d::D1()
  {
    CPoint<1> k;

    ft(-1);
    #pragma omp parallel for private(k)
    for (int i=0; i<m_dim; i++ )
    {
      double tmp = m_out[i][0];
      k = Get_k( i );
      m_out[i][0] = k[0]*m_out[i][1];
      m_out[i][1] = -k[0]*tmp;
    }
    ft(1);
  }

  /**
   * \brief Calculate 2nd derivative of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_1d::D2()
  {
    CPoint<1> k;

    ft(-1);
    #pragma omp parallel for private(k)
    for (int i=0; i<m_dim; i++ )
    {
      k = Get_k( i );
      double f = -k[0]*k[0];
      m_out[i][0] = f*m_out[i][0];
      m_out[i][1] = f*m_out[i][1];
    }
    ft(1);
  }

  /**
   * \brief Performs reordering and rescaling for fourier transformation.
   *
   * @param data Pointer to data to be reordered and rescaled
   * @param d Stepsize
   */
  void cft_1d::fix( fftw_complex *data, const double d )
  {
    double fak = d / sqrt(2.0*M_PI);

    fftw_complex tmp;

    for ( int i=0; i<m_shift_x; i++ )
    {
      memcpy( &tmp, &data[i+m_shift_x], sizeof(fftw_complex) );
      memcpy( &data[i+m_shift_x], &data[i], sizeof(fftw_complex) );
      memcpy( &data[i], &tmp, sizeof(fftw_complex) );
      data[i][0] *= fak;
      data[i][1] *= fak;
      data[i+m_shift_x][0] *= fak;
      data[i+m_shift_x][1] *= fak;
    }

    for ( int i=1; i<m_dim; i+=2 )
    {
      data[i][0] *= -1.0;
      data[i][1] *= -1.0;
    }
  }

  /**
   * \brief Scale data
   *
   * Scale data by factor sx / sqrt(2*Pi)
   *
   * @param data Pointer to data to be scaled
   * @param sx Stepsize
   */
  void cft_1d::scale( fftw_complex *data, const double sx )
  {
    const double fak = sx / sqrt(2.0*M_PI);

    #pragma omp parallel for
    for ( int i=0; i<m_dim; i++ )
    {
      data[i][0] *= fak;
      data[i][1] *= fak;
    }
  }
}
