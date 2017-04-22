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

#include <cstdlib>
#include <cmath>
#include <cstring>
#include "cft_3d.h"

namespace Fourier
{
  /**
   * \brief cft_3d Constructor
   *
   * Complex Fourier Transformation in 3 dimension.
   *
   * @param header Header information to construct cft object
   * @param b Whether inplace transformation is done
   */
  cft_3d::cft_3d( const generic_header &header, bool b, bool f ) : cft_base( header, b, f )
  {
    m_forwardPlan  = fftw_plan_dft_3d( m_dim_x, m_dim_y, m_dim_z, m_in, m_out, FFTW_FORWARD, FFTW_ESTIMATE );
    m_backwardPlan = fftw_plan_dft_3d( m_dim_x, m_dim_y, m_dim_z, m_out, m_in, FFTW_BACKWARD, FFTW_ESTIMATE );

    assert( m_forwardPlan != nullptr );
    assert( m_backwardPlan != nullptr );
  }

  /**
   * \brief Performs Fourier Transformation
   *
   * Performs a Fourier Transformation on members in cft_3d object.
   * Forward FT (isign = -1) transforms data in m_in --> m_out
   * Backward FT (isign = 1) transforms data in m_out --> m_in
   *
   * @param isign Whether forward [isign = -1] or backward [isign = 1]
   fourier transformation is performed.
  */
  void cft_3d::ft( int isign )
  {
    m_isign = isign;
    if ( abs(isign) != 1 ) return;
    if ( isign == -1 )
    {
      fftw_execute( m_forwardPlan );
      if ( m_bfix ) fix( m_out, m_dx, m_dy, m_dz );
      else scale( m_out, m_dx, m_dy, m_dz );
    }
    else
    {
      fftw_execute( m_backwardPlan );
      if ( m_bfix ) fix( m_in, m_dkx, m_dky, m_dkz );
      else scale( m_in, m_dkx, m_dky, m_dkz );
    }
  }

  /**
   * \brief Get k
   *
   * @param[in] i Array Index in x direction
   * @param[in] j Array Index in y direction
   * @param[in] k Array Index in z direction
   * @param[out] t_i Translated Array Index of x dimension
   * @param[out] t_j Translated Array Index of y dimension
   * @param[out] t_k Translated Array Index of z dimension
   * @param[out] k_x x-component of transform variable k
   * @param[out] k_y y-component of transform variable k
   * @param[out] k_z z-component of transform variable k
   */
  void cft_3d::Get_k( const int i, const int j, const int k, int &t_i, int &t_j, int &t_k, double &k_x, double &k_y, double &k_z )
  {
    if ( !m_bfix )
    {
      t_i = (i+m_shift_x)%m_dim_x;
      t_j = (j+m_shift_y)%m_dim_y;
      t_k = (k+m_shift_z)%m_dim_z;
    }
    else
    {
      t_i = i;
      t_j = j;
      t_k = k;
    }

    k_x = m_dkx*double(t_i-m_shift_x);
    k_y = m_dky*double(t_j-m_shift_y);
    k_z = m_dkz*double(t_k-m_shift_z);
  }


  void cft_3d::Get_k( int i, int j, int k, double &k_x, double &k_y, double &k_z )
  {
    if ( !m_bfix )
    {
      i = (i+m_shift_x)%m_dim_x;
      j = (j+m_shift_y)%m_dim_y;
      k = (k+m_shift_z)%m_dim_z;
    }

    k_x = m_dkx*double(i-m_shift_x);
    k_y = m_dky*double(j-m_shift_y);
    k_z = m_dkz*double(k-m_shift_z);
  }
  /**
   * \brief Get x
   *
   * @param l Array Index
   */
  CPoint<3> cft_3d::Get_x( const int64_t l )
  {
    CPoint<3> retval;
    int64_t i = l / m_dim_y / m_dim_z;
    int64_t j = (l - i*m_dim_y*m_dim_z) / m_dim_z;
    int64_t k = l - i*m_dim_y*m_dim_z - j*m_dim_z;
    retval[0] = double(i-m_shift_x)*m_dx;
    retval[1] = double(j-m_shift_y)*m_dy;
    retval[2] = double(k-m_shift_z)*m_dz;
    return retval;
  }

  /**
   * \brief Get k
   *
   * @param l Array Index
   */
  CPoint<3> cft_3d::Get_k( const int64_t l )
  {
    CPoint<3> retval;
    int64_t i = l / m_dim_y / m_dim_z;
    int64_t j = (l - i*m_dim_y*m_dim_z) / m_dim_z;
    int64_t k = l - i*m_dim_y*m_dim_z - j*m_dim_z;
    if ( !m_bfix )
    {
      i = (i+m_shift_x)%m_dim_x;
      j = (j+m_shift_y)%m_dim_y;
      k = (k+m_shift_z)%m_dim_z;
    }

    assert( i < m_dim_x );
    assert( j < m_dim_y );
    assert( k < m_dim_z );

    retval[0] = m_dkx*double(i-m_shift_x);
    retval[1] = m_dky*double(j-m_shift_y);
    retval[2] = m_dkz*double(k-m_shift_z);
    return retval;
  }

  /**
   * \brief Get kx
   *
   * @param[in] i Array Index in x direction
   * @param[out] t_i Translated Array Index of x dimension
   * @param[out] k_x x-component of transform variable k
   */
  void cft_3d::Get_kx( const int i, int &t_i, double &k_x )
  {
    if ( !m_bfix )
      t_i = (i+m_shift_x)%m_dim_x;
    else
      t_i = i;

    k_x = m_dkx*double(t_i-m_shift_x);
  }

  /**
   * \brief Get ky
   *
   * @param[in] j Array Index in x direction
   * @param[out] t_j Translated Array Index of x dimension
   * @param[out] k_y y-component of transform variable k
   */
  void cft_3d::Get_ky( const int j, int &t_j, double &k_y )
  {
    if ( !m_bfix )
      t_j = (j+m_shift_y)%m_dim_y;
    else
      t_j = j;

    k_y = m_dky*double(t_j-m_shift_y);
  }

  /**
   * \brief Get kz
   *
   * @param[in] k Array Index in z direction
   * @param[out] t_k Translated Array Index of z dimension
   * @param[out] k_z z-component of transform variable k
   */
  void cft_3d::Get_kz( const int k, int &t_k, double &k_z )
  {
    if ( !m_bfix )
      t_k = (k+m_shift_z)%m_dim_z;
    else
      t_k = k;

    k_z = m_dkz*double(t_k-m_shift_z);
  }

  /**
   * \brief Get kx
   *
   * @param[in] i Array Index in x direction
   * @return x-component of transform variable k
   */
  double cft_3d::Get_kx( const int i )
  {
    int t_i;
    if ( !m_bfix )
      t_i = (i+m_shift_x)%m_dim_x;
    else
      t_i = i;
    return m_dkx*double(t_i-m_shift_x);
  }

  /**
   * \brief Get ky
   *
   * @param[in] j Array Index in y direction
   * @return y-component of transform variable k
   */
  double cft_3d::Get_ky( const int j )
  {
    int t_j;
    if ( !m_bfix )
      t_j = (j+m_shift_y)%m_dim_y;
    else
      t_j = j;
    return m_dky*double(t_j-m_shift_y);
  }

  /**
   * \brief Get kz
   *
   * @param[in] k Array Index in z direction
   * @return z-component of transform variable k
   */
  double cft_3d::Get_kz( const int k )
  {
    int t_k;
    if ( !m_bfix )
      t_k = (k+m_shift_z)%m_dim_z;
    else
      t_k = k;
    return m_dkz*double(t_k-m_shift_z);
  }

  /**
   * \brief Calculate 1st derivative with respect to x of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_3d::Diff_x()
  {
    int ijk, i2;
    double kx, tmp1;

    bool b_oldval = m_bfix;
    m_bfix = false;

    this->ft(-1);
    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk,i2,kx,tmp1)
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_dim_z; k++ )
        {
          this->Get_kx( i, i2, kx );
          ijk = k+m_dim_z*(j+m_dim_y*i);

          tmp1 = m_out[ijk][0];
          m_out[ijk][0] = -kx*m_out[ijk][1];
          m_out[ijk][1] = kx*tmp1;
        }
      }
    }
    this->ft(1);
    m_bfix = b_oldval;
  }

  /**
   * \brief Calculate 1st derivative with respect to y of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_3d::Diff_y()
  {
    int ijk, j2;
    double ky, tmp1;

    bool b_oldval = m_bfix;
    m_bfix = false;

    this->ft(-1);
    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk,j2,ky,tmp1)
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_dim_z; k++ )
        {
          this->Get_ky( j, j2, ky );
          ijk = k+m_dim_z*(j+m_dim_y*i);

          tmp1 = m_out[ijk][0];
          m_out[ijk][0] = -ky*m_out[ijk][1];
          m_out[ijk][1] = ky*tmp1;
        }
      }
    }
    this->ft(1);
    m_bfix = b_oldval;
  }

  /**
   * \brief Calculate 1st derivative with respect to z of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_3d::Diff_z()
  {
    int ijk, k2;
    double kz, tmp1;

    bool b_oldval = m_bfix;
    m_bfix = false;

    this->ft(-1);
    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk,k2,kz,tmp1)
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_dim_z; k++ )
        {
          this->Get_kz( k, k2, kz );
          ijk = k+m_dim_z*(j+m_dim_y*i);

          tmp1 = m_out[ijk][0];
          m_out[ijk][0] = -kz*m_out[ijk][1];
          m_out[ijk][1] = kz*tmp1;
        }
      }
    }
    this->ft(1);
    m_bfix = b_oldval;
  }

  /**
   * \brief Calculate 2nd derivative with respect to x of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_3d::Diff_xx()
  {
    int ijk, i2;
    double kx;

    bool b_oldval = m_bfix;
    m_bfix = false;

    this->ft(-1);
    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk,i2,kx)
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_dim_z; k++ )
        {
          this->Get_kx( i, i2, kx );
          kx = -kx*kx;
          ijk = k+m_dim_z*(j+m_dim_y*i);

          m_out[ijk][0] *= kx;
          m_out[ijk][1] *= kx;
        }
      }
    }
    this->ft(1);
    m_bfix = b_oldval;
  }

  /**
   * \brief Calculate 2nd derivative with respect to y of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_3d::Diff_yy()
  {
    int ijk, j2;
    double ky;

    bool b_oldval = m_bfix;
    m_bfix = false;

    this->ft(-1);
    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk,j2,ky)
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_dim_z; k++ )
        {
          this->Get_ky( j, j2, ky );
          ky = -ky*ky;
          ijk = k+m_dim_z*(j+m_dim_y*i);

          m_out[ijk][0] *= ky;
          m_out[ijk][1] *= ky;
        }
      }
    }
    this->ft(1);
    m_bfix = b_oldval;
  }

  /**
   * \brief Calculate 2nd derivative with respect to z of data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_3d::Diff_zz()
  {
    int ijk, k2;
    double kz;

    bool b_oldval = m_bfix;
    m_bfix = false;

    this->ft(-1);
    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk,k2,kz)
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_dim_z; k++ )
        {
          this->Get_kz( k, k2, kz );
          kz = -kz*kz;
          ijk = k+m_dim_z*(j+m_dim_y*i);

          m_out[ijk][0] *= kz;
          m_out[ijk][1] *= kz;
        }
      }
    }
    this->ft(1);
    m_bfix = b_oldval;
  }

  /**
   * \brief Applies Laplace Operator on data in m_in
   *
   *  Differentiation is done via fourier transformation method.
   */
  void cft_3d::Laplace()
  {
    int ijk, i2, j2, k2;
    double kx, ky, kz, tmp1;

    bool b_oldval = m_bfix;
    m_bfix = false;

    this->ft(-1);
    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(i2,j2,k2,kx,ky,kz,tmp1)
    for ( int i=0; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_dim_z; k++ )
        {
          this->Get_k( i, j, k, i2, j2, k2, kx, ky, kz );
          ijk = k+m_dim_z*(j+m_dim_y*i);

          tmp1 = -kx*kx-ky*ky-kz*kz;

          m_out[ijk][0] *= tmp1;
          m_out[ijk][1] *= tmp1;
        }
      }
    }
    this->ft(1);
    m_bfix = b_oldval;
  }

  /**
   * \brief Performs reordering and rescaling after fourier transformation.
   *
   * @param data Pointer to data to be reordered and rescaled
   * @param sx Stepsize in x direction
   * @param sy Stepsize in y direction
   * @param sz Stepsize in z direction
   */
  void cft_3d::fix( fftw_complex *data, const double sx, const double sy, const double sz )
  {
    int ijk_1, ijk_2;
    double fak2, tmp;
    const double fak = sx * sy * sz / pow(2*M_PI,1.5);

    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk_1,ijk_2,fak2,tmp)
    for ( int i=0; i<m_shift_x; i++ )
    {
      for ( int j=0; j<m_shift_y; j++ )
      {
        for ( int k=0; k<m_shift_z; k++ )
        {
          ijk_1 = k+m_dim_z*(j+m_dim_y*i);
          ijk_2 = (k+m_shift_z)+m_dim_z*((j+m_shift_y)+m_dim_y*(i+m_shift_x));

          if ( (i+j+k)%2 == 1 ) fak2 = -fak;
          else fak2 = fak;

          tmp            = data[ijk_1][0];
          data[ijk_1][0] = data[ijk_2][0] * fak2;
          data[ijk_2][0] = tmp * fak2;

          tmp            = data[ijk_1][1];
          data[ijk_1][1] = data[ijk_2][1] * fak2;
          data[ijk_2][1] = tmp * fak2;
        }
      }
    }

    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk_1,ijk_2,fak2,tmp)
    for ( int i=m_shift_x; i<m_dim_x; i++ )
    {
      for ( int j=0; j<m_shift_y; j++ )
      {
        for ( int k=0; k<m_shift_z; k++ )
        {
          ijk_1 = k+m_dim_z*(j+m_dim_y*i);
          ijk_2 = (k+m_shift_z)+m_dim_z*((j+m_shift_y)+m_dim_y*(i-m_shift_x));

          if ( (i+j+k)%2 == 1 ) fak2 = -fak;
          else fak2 = fak;

          tmp            = data[ijk_1][0];
          data[ijk_1][0] = data[ijk_2][0] * fak2;
          data[ijk_2][0] = tmp * fak2;

          tmp            = data[ijk_1][1];
          data[ijk_1][1] = data[ijk_2][1] * fak2;
          data[ijk_2][1] = tmp * fak2;
        }
      }
    }

    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk_1,ijk_2,fak2,tmp)
    for ( int i=0; i<m_shift_x; i++ )
    {
      for ( int j=m_shift_y; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_shift_z; k++ )
        {
          ijk_1 = k+m_dim_z*(j+m_dim_y*i);
          ijk_2 = (k+m_shift_z)+m_dim_z*((j-m_shift_y)+m_dim_y*(i+m_shift_x));

          if ( (i+j+k)%2 == 1 ) fak2 = -fak;
          else fak2 = fak;

          tmp            = data[ijk_1][0];
          data[ijk_1][0] = data[ijk_2][0] * fak2;
          data[ijk_2][0] = tmp * fak2;

          tmp            = data[ijk_1][1];
          data[ijk_1][1] = data[ijk_2][1] * fak2;
          data[ijk_2][1] = tmp * fak2;
        }
      }
    }

    #pragma omp parallel for schedule(dynamic,1) collapse(3) private(ijk_1,ijk_2,fak2,tmp)
    for ( int i=m_shift_x; i<m_dim_x; i++ )
    {
      for ( int j=m_shift_y; j<m_dim_y; j++ )
      {
        for ( int k=0; k<m_shift_z; k++ )
        {
          ijk_1 = k+m_dim_z*(j+m_dim_y*i);
          ijk_2 = (k+m_shift_z)+m_dim_z*((j-m_shift_y)+m_dim_y*(i-m_shift_x));

          if ( (i+j+k)%2 == 1 ) fak2 = -fak;
          else fak2 = fak;

          tmp            = data[ijk_1][0];
          data[ijk_1][0] = data[ijk_2][0] * fak2;
          data[ijk_2][0] = tmp * fak2;

          tmp            = data[ijk_1][1];
          data[ijk_1][1] = data[ijk_2][1] * fak2;
          data[ijk_2][1] = tmp * fak2;
        }
      }
    }
  }

  /**
   * \brief Scale data
   *
   * Scale data by factor sx*sy / sqrt(2*Pi)
   *
   * @param data Pointer to data to be scaled
   * @param sx Stepsize in x direction
   * @param sy Stepsize in y direction
   * @param sz Stepsize in z direction
   */
  void cft_3d::scale( fftw_complex *data, const double sx, const double sy, const double sz )
  {
    const double fak = sx * sy * sz / pow(2*M_PI,1.5);

    #pragma omp parallel for
    for ( int64_t i=0; i<m_dim; i++ )
    {
      data[i][0] *= fak;
      data[i][1] *= fak;
    }
  }
}
