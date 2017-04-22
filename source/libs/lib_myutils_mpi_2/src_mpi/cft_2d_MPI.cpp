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
#include "cft_2d_MPI.h"

using namespace std;

namespace MPI
{
  namespace Fourier
  {
    /**
     * \brief cft_2d_MPI Constructor
     *
     * Complex Fourier Transformation in 2 dimension. MPI Version.
     *
     * @param header Header information to construct cft object
     */
    cft_2d_MPI::cft_2d_MPI( generic_header *header ) : cft_base_MPI<2>(header)
    {

      ptrdiff_t alloc_local = fftw_mpi_local_size_2d_transposed( m_dimX, m_dimY, MPI_COMM_WORLD, &m_loc_dimX, &m_loc_start_dimX, &m_loc_dimY, &m_loc_start_dimY );
      //printf("Rank %d || dimX: %ld | dimY: %ld | locX: %ld | loc_start_X: %ld | locY: %ld | loc_start_Y: %ld\n", m_rank, m_dimX, m_dimY, m_loc_dimX,m_loc_start_dimX,m_loc_dimY,m_loc_start_dimY);

      m_data = fftw_alloc_complex(alloc_local);

      m_forwardPlan = fftw_mpi_plan_dft_2d( m_dimX, m_dimY, m_data, m_data, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_OUT);
      m_backwardPlan = fftw_mpi_plan_dft_2d( m_dimX, m_dimY, m_data, m_data, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_MPI_TRANSPOSED_IN);

      m_offset_rs = m_loc_start_dimX*m_dimY;
      m_offset_fs = m_loc_start_dimY*m_dimX;
      m_loc_n_rs = m_loc_dimX*m_dimY;
      m_loc_n_fs = m_loc_dimY*m_dimX;
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
    void cft_2d_MPI::ft( const int isign )
    {
      double fak;

      if ( isign == -1 )
      {
        fftw_execute( m_forwardPlan );

        m_fs = true;
        fak = 0.5*m_dx*m_dy/M_PI;

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
        fak = 0.5*m_dkx*m_dky/M_PI;

        for ( ptrdiff_t l=0; l<m_loc_n_rs; l++ )
        {
          m_data[l][0] *= fak;
          m_data[l][1] *= fak;
        }
      }
      MPI_Barrier( MPI_COMM_WORLD );
    }

    /**
     * \brief Get x value
     *
     * @param loc Local linear Array Index
     */
    CPoint<2> cft_2d_MPI::Get_x( const ptrdiff_t loc )
    {
      CPoint<2> retval;
      ptrdiff_t i = (m_offset_rs+loc) / m_dimY;
      ptrdiff_t j = (m_offset_rs+loc) - i*m_dimY;

      retval[0] = double(i-m_shift_x)*m_dx;
      retval[1] = double(j-m_shift_y)*m_dy;
      return retval;
    }

    /**
     * \brief Get k value
     *
     * @param loc Local linear Array Index
     */
    CPoint<2> cft_2d_MPI::Get_k( const ptrdiff_t loc )
    {
      CPoint<2> retval;
      ptrdiff_t i = (m_offset_fs+loc) / m_dimX;
      ptrdiff_t j = (m_offset_fs+loc) - i*m_dimX;

      retval[0] = m_dkx*double(((j+m_shift_x)%m_dimX)-m_shift_x);
      retval[1] = m_dky*double(((i+m_shift_y)%m_dimY)-m_shift_y);
      return retval;
    }

    /**
     * \brief Get k
     *
     * @param[in] i Array Index in x direction
     * @param[in] j Array Index in y direction
     * @param[out] t_i Translated Array Index of x dimension
     * @param[out] t_j Translated Array Index of y dimension
     * @param[out] k_x x-component of transform variable k
     * @param[out] k_y y-component of transform variable k
     */
    void cft_2d_MPI::Get_k( const ptrdiff_t i, const ptrdiff_t j, ptrdiff_t &t_i, ptrdiff_t &t_j, double &k_x, double &k_y )
    {
      t_i = (i+m_shift_x)%m_dimX;
      k_x = m_dkx*double(t_i-m_shift_x);

      t_j = (j+m_loc_start_dimY+m_shift_y)%m_dimY;
      k_y = m_dky*double(t_j-m_shift_y);
    }

    /**
     * \brief Get array indices of global matrix in real space
     *
     * @param loc Local linear Array Index
     * @param i Global Array Index of x dimension
     * @param j Global Array Index of y dimension
     */
    void cft_2d_MPI::local_to_global_rs( const ptrdiff_t loc, ptrdiff_t &i, ptrdiff_t &j )
    {
      i = (m_offset_rs+loc) / m_dimY;
      j = (m_offset_rs+loc) - i*m_dimY;
    }

    void cft_2d_MPI::local_to_global_rs( const ptrdiff_t loc, std::vector<ptrdiff_t> &retval )
    {
      retval[0] = (m_offset_rs+loc) / m_dimY;
      retval[1] = (m_offset_rs+loc) - retval[0]*m_dimY;
    }

    /**
     * \brief Get array indices of global matrix in fourier space
     *
     * @param loc Local linear Array Index
     * @param i Global Array Index of x dimension
     * @param j Global Array Index of y dimension
     */
    void cft_2d_MPI::local_to_global_fs( const ptrdiff_t loc, ptrdiff_t &i, ptrdiff_t &j )
    {
      i = (m_offset_fs+loc) / m_dimX;
      j = (m_offset_fs+loc) - i*m_dimX;
    }

    void cft_2d_MPI::local_to_global_fs( const ptrdiff_t loc, std::vector<ptrdiff_t> &retval )
    {
      retval[0] = (m_offset_fs+loc) / m_dimX;
      retval[1] = (m_offset_fs+loc) - retval[0]*m_dimX;
    }

    /**
     * \brief Get x value from local linear array index
     *
     * @param loc Local linear Array Index
     */
    void cft_2d_MPI::local_to_x( const ptrdiff_t loc, CPoint<2> &x )
    {
      ptrdiff_t i = (m_offset_rs+loc) / m_dimY;
      ptrdiff_t j = (m_offset_rs+loc) - i*m_dimY;

      x[0] = double(i-m_shift_x)*m_dx;
      x[1] = double(j-m_shift_y)*m_dy;
    }

    /**
     * \brief Get k value from local linear array index
     *
     * @param loc Local linear Array Index
     */
    void cft_2d_MPI::local_to_k( const ptrdiff_t loc, CPoint<2> &k )
    {
      ptrdiff_t i = (m_offset_fs+loc) / m_dimX;
      ptrdiff_t j = (m_offset_fs+loc) - i*m_dimX;

      k[0] = m_dkx*double(((j+m_shift_x)%m_dimX)-m_shift_x);
      k[1] = m_dky*double(((i+m_shift_y)%m_dimY)-m_shift_y);
    }

    /**
     * \brief Applies Laplace Operator on data in m_in
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_2d_MPI::Laplace()
    {
      ft(-1);

      CPoint<2> k;
      for ( ptrdiff_t i=0; i<Get_loc_n_fs(); i++ )
      {
        local_to_k( i, k );

        double fak = -(k*k);
        m_data[i][0] *= fak;
        m_data[i][1] *= fak;
      }

      ft(1);
    }

    /**
     * \brief Calculate 1st derivative with respect to x of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_2d_MPI::D_x()
    {
      ft(-1);

      CPoint<2> k;
      for ( ptrdiff_t i=0; i<Get_loc_n_fs(); i++ )
      {
        local_to_k( i, k );

        double tmp1 = m_data[i][0];
        m_data[i][0] = -k[0]*m_data[i][1];
        m_data[i][1] = k[0]*tmp1;
      }

      ft(1);
    }

    /**
     * \brief Calculate 1st derivative with respect to y of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_2d_MPI::D_y()
    {
      ft(-1);

      CPoint<2> k;
      for ( ptrdiff_t i=0; i<Get_loc_n_fs(); i++ )
      {
        local_to_k( i, k );

        double tmp1 = m_data[i][0];
        m_data[i][0] = -k[1]*m_data[i][1];
        m_data[i][1] = k[1]*tmp1;
      }

      ft(1);
    }

    /**
     * \brief Calculate 2nd derivative with respect to x of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_2d_MPI::D_xx()
    {
      ft(-1);

      CPoint<2> k;
      for ( ptrdiff_t i=0; i<Get_loc_n_fs(); i++ )
      {
        local_to_k( i, k );

        double fak = -k[0]*k[0];
        m_data[i][0] *= fak;
        m_data[i][1] *= fak;
      }

      ft(1);
    }

    /**
     * \brief Calculate 2nd derivative with respect to y of data
     *
     *  Differentiation is done via fourier transformation method.
     */
    void cft_2d_MPI::D_yy()
    {
      ft(-1);

      CPoint<2> k;
      for ( ptrdiff_t i=0; i<Get_loc_n_fs(); i++ )
      {
        local_to_k( i, k );

        double fak = -k[1]*k[1];
        m_data[i][0] *= fak;
        m_data[i][1] *= fak;
      }

      ft(1);
    }
  }
}
