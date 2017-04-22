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
#include "fftw3-mpi.h"
#include "cft_3d_MPI.h"

using namespace std;

namespace MPI
{
  namespace Fourier
  {
    /**
     * \brief cft_3d Constructor
     *
     * Complex Fourier Transformation in 3 dimension. MPI Version
     *
     * @param header Header information to construct cft object
     */
    cft_3d_MPI::cft_3d_MPI( generic_header *header ) : cft_base_MPI<3>(header)
    {
      ptrdiff_t alloc_local = fftw_mpi_local_size_3d_transposed( m_dimX, m_dimY, m_dimZ, MPI_COMM_WORLD, &m_loc_dimX, &m_loc_start_dimX, &m_loc_dimY, &m_loc_start_dimY );

      m_data = fftw_alloc_complex(alloc_local);

      m_forwardPlan = fftw_mpi_plan_dft_3d( m_dimX, m_dimY, m_dimZ, m_data, m_data, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_OUT);
      m_backwardPlan = fftw_mpi_plan_dft_3d( m_dimX, m_dimY, m_dimZ, m_data, m_data, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);

      m_offset_rs = m_loc_start_dimX*m_dimY*m_dimZ;
      m_offset_fs = m_loc_start_dimY*m_dimX*m_dimZ;
      m_loc_n_rs = m_loc_dimX*m_dimY*m_dimZ;
      m_loc_n_fs = m_loc_dimY*m_dimX*m_dimZ;
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
    void cft_3d_MPI::ft( const int isign )
    {
      double fak;

      if ( isign == -1 )
      {
        fftw_execute( m_forwardPlan );

        m_fs = true;
        fak = m_dx*m_dy*m_dz/pow( 2.0*M_PI, 1.5 );

        for ( ptrdiff_t l=0; l<m_loc_n_fs; l++ )
        {
          m_data[l][0] *= fak;
          m_data[l][1] *= fak;
        }
      }
      else
      {
        fftw_execute( m_backwardPlan );

        m_fs = false;
        fak = m_dkx*m_dky*m_dkz/pow( 2.0*M_PI, 1.5 );

        for ( ptrdiff_t l=0; l<m_loc_n_rs; l++ )
        {
          m_data[l][0] *= fak;
          m_data[l][1] *= fak;
        }
      }
    }

    /**
     * \brief Get x
     *
     * @param loc Local Linear Array Index
     */
    CPoint<3> cft_3d_MPI::Get_x( const ptrdiff_t loc )
    {
      CPoint<3> retval;
      ptrdiff_t L = m_offset_rs+loc;
      ptrdiff_t i = L / m_dimYZ;
      ptrdiff_t j = L / m_dimZ - i*m_dimY;
      ptrdiff_t k = L - i*m_dimYZ - j*m_dimZ;

      retval[0] = double(i-m_shift_x)*m_dx;
      retval[1] = double(j-m_shift_y)*m_dy;
      retval[2] = double(k-m_shift_z)*m_dz;
      return retval;
    }

    /**
     * \brief Get k
     *
     * @param loc Local Linear Array Index
     */
    CPoint<3> cft_3d_MPI::Get_k( const ptrdiff_t loc )
    {
      CPoint<3> retval;
      ptrdiff_t L = m_offset_fs+loc;
      ptrdiff_t i = L / m_dimXZ;
      ptrdiff_t j = L / m_dimZ - i*m_dimX;
      ptrdiff_t k = L - i*m_dimXZ - j*m_dimZ;

      retval[0] = m_dkx*double(((j+m_shift_x)%m_dimX)-m_shift_x);
      retval[1] = m_dky*double(((i+m_shift_y)%m_dimY)-m_shift_y);
      retval[2] = m_dkz*double(((k+m_shift_z)%m_dimZ)-m_shift_z);
      return retval;
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
    void cft_3d_MPI::Get_k( const ptrdiff_t i, const ptrdiff_t j, const ptrdiff_t k, ptrdiff_t &t_i, ptrdiff_t &t_j, ptrdiff_t &t_k, double &k_x, double &k_y, double &k_z )
    {
      t_i = (i+m_shift_x)%m_dimX;
      k_x = m_dkx*double(t_i-m_shift_x);

      t_j = (j+m_loc_start_dimY+m_shift_y)%m_dimY;
      k_y = m_dky*double(t_j-m_shift_y);

      t_k = (k+m_shift_z)%m_dimZ;
      k_z = m_dkz*double(t_k-m_shift_z);
    }

    /**
     * \brief Get kz
     *
     * @param[in] k Array Index in z direction
     * @param[out] t_k Translated Array Index of z dimension
     * @param[out] k_z z-component of transform variable k
     */
    void cft_3d_MPI::Get_kz( const ptrdiff_t k, ptrdiff_t &t_k, double &k_z )
    {
      t_k = (k+m_shift_z)%m_dimZ;
      k_z = m_dkz*double(t_k-m_shift_z);
    }

    /**
     * \brief Get kz
     *
     * @param[in] k Array Index in z direction
     * @return z-component of transform variable k
     */
    double cft_3d_MPI::Get_kz( const ptrdiff_t k )
    {
      return m_dkz*double(((k+m_shift_z)%m_dimZ)-m_shift_z);
    }

    /**
     * \brief Get z
     *
     * @param[in] k Array Index in z direction
     * @return z-coordinate
     */
    double cft_3d_MPI::Get_z( const ptrdiff_t k )
    {
      return double(k-m_shift_z)*m_dz;
    }

    /**
     * \brief Get array indices of global matrix in real space
     *
     * @param loc Local linear Array Index
     * @param i Global Array Index of x dimension
     * @param j Global Array Index of y dimension
     */
    void cft_3d_MPI::local_to_global_rs( const ptrdiff_t loc, ptrdiff_t &i, ptrdiff_t &j, ptrdiff_t &k )
    {
      ptrdiff_t L = m_offset_rs+loc;
      i = L / m_dimYZ;
      j = L / m_dimZ - i*m_dimY;
      k = L - i*m_dimYZ - j*m_dimZ;
    }

    void cft_3d_MPI::local_to_global_rs( const ptrdiff_t loc, std::vector<ptrdiff_t> &retval )
    {
      ptrdiff_t L = m_offset_rs+loc;
      retval[0] = L / m_dimYZ;
      retval[1] = L / m_dimZ -  retval[0]*m_dimY;
      retval[2] = L -  retval[0]*m_dimYZ -  retval[1]*m_dimZ;
    }

    /**
     * \brief Get array indices of global matrix in fourier space
     *
     * @param loc Local linear Array Index
     * @param i Global Array Index of x dimension
     * @param j Global Array Index of y dimension
     */
    void cft_3d_MPI::local_to_global_fs( const ptrdiff_t loc, ptrdiff_t &i, ptrdiff_t &j, ptrdiff_t &k )
    {
      ptrdiff_t L = m_offset_fs+loc;
      i = L / m_dimXZ;
      j = L / m_dimZ - i*m_dimX;
      k = L - i*m_dimXZ - j*m_dimZ;
    }

    void cft_3d_MPI::local_to_global_fs( const ptrdiff_t loc, std::vector<ptrdiff_t> &retval )
    {
      ptrdiff_t L = m_offset_fs+loc;
      retval[0] = L / m_dimXZ;
      retval[1] = L / m_dimZ - retval[0]*m_dimX;
      retval[2] = L - retval[0]*m_dimXZ - retval[1]*m_dimZ;
    }

    /**
     * \brief Get x value from local linear array index
     *
     * @param loc Local linear Array Index
     */
    void cft_3d_MPI::local_to_x( const ptrdiff_t loc, CPoint<3> &x )
    {
      ptrdiff_t L = m_offset_rs+loc;
      ptrdiff_t i = L / m_dimYZ;
      ptrdiff_t j = L / m_dimZ - i*m_dimY;
      ptrdiff_t k = L - i*m_dimYZ - j*m_dimZ;

      x[0] = double(i-m_shift_x)*m_dx;
      x[1] = double(j-m_shift_y)*m_dy;
      x[2] = double(k-m_shift_z)*m_dz;
    }

    /**
     * \brief Get k value from local linear array index
     *
     * @param loc Local linear Array Index
     */
    void cft_3d_MPI::local_to_k( const ptrdiff_t loc, CPoint<3> &x )
    {
      ptrdiff_t L = m_offset_fs+loc;
      ptrdiff_t i = L / m_dimXZ;
      ptrdiff_t j = L / m_dimZ - i*m_dimX;
      ptrdiff_t k = L - i*m_dimXZ - j*m_dimZ;

      x[0] = m_dkx*double(((j+m_shift_x)%m_dimX)-m_shift_x);
      x[1] = m_dky*double(((i+m_shift_y)%m_dimY)-m_shift_y);
      x[2] = m_dkz*double(((k+m_shift_z)%m_dimZ)-m_shift_z);
    }

    /**
     * \brief Applies Laplace Operator on data in m_in
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_3d_MPI::Laplace()
    {
      ft(-1);

      ptrdiff_t i, j, k, ijk, ti, tj, tk;
      double kx, ky, kz, fak;

      for ( j=0; j<m_loc_dimY; j++ )
      {
        Get_ky(j, tj, ky);
        for ( i=0; i<m_dimX; i++ )
        {
          Get_kx(i, ti, kx);

          for ( k=0; k<m_dimZ; k++ )
          {
            Get_kz(k, tk, kz);

            ijk = k+m_dimZ*(i+m_dimX*j);

            fak = -kx*kx-ky*ky-kz*kz;
            m_data[ijk][0] *= fak;
            m_data[ijk][1] *= fak;
          }
        }
      }

      ft(1);
    }

    /**
     * \brief Calculate 1st derivative with respect to x of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_3d_MPI::D_x()
    {
      ft(-1);

      ptrdiff_t i, j, k, ijk, ti;
      double kx, tmp1;

      for ( j=0; j<m_loc_dimY; j++ )
      {
        for ( i=0; i<m_dimX; i++ )
        {
          Get_kx(i, ti, kx);
          for ( k=0; k<m_dimZ; k++ )
          {
            ijk = k+m_dimZ*(i+m_dimX*j);

            tmp1 = m_data[ijk][0];
            m_data[ijk][0] = -kx*m_data[ijk][1];
            m_data[ijk][1] = kx*tmp1;
          }
        }
      }

      ft(1);
    }

    /**
     * \brief Calculate 1st derivative with respect to y of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_3d_MPI::D_y()
    {
      ft(-1);

      ptrdiff_t i, j, k, ijk, tj;
      double ky, tmp1;

      for ( j=0; j<m_loc_dimY; j++ )
      {
        Get_ky(j, tj, ky);
        for ( i=0; i<m_dimX; i++ )
        {
          for ( k=0; k<m_dimZ; k++ )
          {
            ijk = k+m_dimZ*(i+m_dimX*j);

            tmp1 = m_data[ijk][0];
            m_data[ijk][0] = -ky*m_data[ijk][1];
            m_data[ijk][1] = ky*tmp1;
          }
        }
      }

      ft(1);
    }

    /**
     * \brief Calculate 1st derivative with respect to z of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_3d_MPI::D_z()
    {
      ft(-1);

      ptrdiff_t i, j, k, ijk, tk;
      double kz, tmp1;

      for ( j=0; j<m_loc_dimY; j++ )
      {
        for ( i=0; i<m_dimX; i++ )
        {
          for ( k=0; k<m_dimZ; k++ )
          {
            Get_kz(k, tk, kz);

            ijk = k+m_dimZ*(i+m_dimX*j);

            tmp1 = m_data[ijk][0];
            m_data[ijk][0] = -kz*m_data[ijk][1];
            m_data[ijk][1] = kz*tmp1;
          }
        }
      }

      ft(1);
    }

    /**
     * \brief Calculate 2nd derivative with respect to x of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_3d_MPI::D_xx()
    {
      ft(-1);

      ptrdiff_t i, j, k, ijk, ti;
      double kx, fak;

      for ( j=0; j<m_loc_dimY; j++ )
      {
        for ( i=0; i<m_dimX; i++ )
        {
          Get_kx(i, ti, kx);

          for ( k=0; k<m_dimZ; k++ )
          {
            ijk = k+m_dimZ*(i+m_dimX*j);

            fak = -kx*kx;
            m_data[ijk][0] *= fak;
            m_data[ijk][1] *= fak;
          }
        }
      }

      ft(1);
    }

    /**
     * \brief Calculate 2nd derivative with respect to y of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_3d_MPI::D_yy()
    {
      ft(-1);

      ptrdiff_t i, j, k, ijk, tj;
      double ky, fak;

      for ( j=0; j<m_loc_dimY; j++ )
      {
        Get_ky(j, tj, ky);
        for ( i=0; i<m_dimX; i++ )
        {
          for ( k=0; k<m_dimZ; k++ )
          {
            ijk = k+m_dimZ*(i+m_dimX*j);

            fak = -ky*ky;
            m_data[ijk][0] *= fak;
            m_data[ijk][1] *= fak;
          }
        }
      }

      ft(1);
    }

    /**
     * \brief Calculate 2nd derivative with respect to z of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_3d_MPI::D_zz()
    {
      ft(-1);

      ptrdiff_t i, j, k, ijk, tk;
      double kz, fak;

      for ( j=0; j<m_loc_dimY; j++ )
      {
        for ( i=0; i<m_dimX; i++ )
        {
          for ( k=0; k<m_dimZ; k++ )
          {
            Get_kz(k, tk, kz);

            ijk = k+m_dimZ*(i+m_dimX*j);

            fak = -kz*kz;
            m_data[ijk][0] *= fak;
            m_data[ijk][1] *= fak;
          }
        }
      }

      ft(1);
    }
  }
}

