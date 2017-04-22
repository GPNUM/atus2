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
#include "fftw3.h"
#include <cmath>
#include "CPoint.h"
#include "my_structs.h"
#include "CPoint.h"

#pragma once

namespace MPI { namespace Fourier
{
  template <int dim>
  class cft_base_MPI
  {
  public:
    /**
    * \brief cft_base_mpi Constructor
    */
    cft_base_MPI( const generic_header* header ) : 
      m_loc_dimX(0), 
      m_loc_dimY(0),
      m_loc_start_dimX(0),
      m_loc_start_dimY(0),
      m_offset_rs(0),
      m_offset_fs(0),
      m_loc_n_rs(0),
      m_loc_n_fs(0)
    {
      m_header = *header;
      m_time = &header->t;

      m_dimX = m_header.nDimX;
      m_dimY = m_header.nDimY;
      m_dimZ = m_header.nDimZ;

      m_dx = m_header.dx;
      m_dy = m_header.dy;
      m_dz = m_header.dz;

      m_dkx = m_header.dkx;
      m_dky = m_header.dky;
      m_dkz = m_header.dkz;

      if( m_dimY==0 ) m_dimY=1;
      if( m_dimZ==0 ) m_dimZ=1;
      m_dimXZ = m_dimX*m_dimZ;
      m_dimYZ = m_dimY*m_dimZ;
      m_fs = false;
      m_shift_x = m_dimX / 2;
      m_shift_y = m_dimY / 2;
      m_shift_z = m_dimZ / 2;

      MPI_Comm_size( MPI_COMM_WORLD, &m_nprocs );
      MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
    }

    /**
    *  \brief cft_base_mpi Destructor
    */
    virtual ~cft_base_MPI()
    {
      fftw_destroy_plan( m_forwardPlan );
      fftw_destroy_plan( m_backwardPlan );
      fftw_free( m_data );
    }

    virtual void ft(int)=0;
    virtual CPoint<dim> Get_k(const ptrdiff_t)=0;
    virtual CPoint<dim> Get_x(const ptrdiff_t)=0;

    fftw_complex * Get_p2_Data() { return m_data; };

    /**
    * \brief Get rank of the calling process in the group of comm
    */
    int Get_Rank()   { return m_rank; };

    /**
    * \brief Get number of processes in the group of comm
    */
    int Get_nprocs() { return m_nprocs; };

    /**
    * \brief Get local number of sampling points in x direction
    */
    ptrdiff_t Get_loc_dimX() { return m_loc_dimX; };

    /**
    * \brief Get local number of sampling points in y direction
    */
    ptrdiff_t Get_loc_dimY() { return m_loc_dimY; };

    /**
    * \brief Get index offset of local x-dimensional matrix part
    */
    ptrdiff_t Get_loc_start_dimX() { return m_loc_start_dimX; };

    /**
    * \brief Get index offset of local y-dimensional matrix part
    */
    ptrdiff_t Get_loc_start_dimY() { return m_loc_start_dimY; };

    /**
    * \brief Get Global number of sampling points in x-direction
    */
    ptrdiff_t Get_dimX() { return m_dimX; };

    /**
    * \brief Get Global number of sampling points in y-direction
    */
    ptrdiff_t Get_dimY() { return m_dimY; };

    /**
    * \brief Get Global number of sampling points in z-direction
    */
    ptrdiff_t Get_dimZ() { return m_dimZ; };

    /**
    * \brief Get offset in real-space
    */
    ptrdiff_t Get_offset_rs() { return m_offset_rs; };

    /**
    * \brief Get offset in fourier-space
    */
    ptrdiff_t Get_offset_fs() { return m_offset_fs; };

    /**
    * \brief Get local number of elements in real space
    */
    ptrdiff_t Get_loc_n_rs() { return m_loc_n_rs; };

    /**
    * \brief Get local number of elements in fourier space
    */
    ptrdiff_t Get_loc_n_fs() { return m_loc_n_fs; };

    /**
    * \brief Get stepsize in x-direction
    */
    double Get_dx() { return m_header.dx; };

    /**
    * \brief Get stepsize in y-direction
    */
    double Get_dy() { return m_header.dy; };

    /**
    * \brief Get stepsize in z-direction
    */
    double Get_dz() { return m_header.dz; };

    /**
    * \brief Get stepsize in kx-direction
    */
    double Get_dkx() { return m_header.dkx; };

    /**
    * \brief Get stepsize in ky-direction
    */
    double Get_dky() { return m_header.dky; };

    /**
    * \brief Get stepsize in kz-direction
    */
    double Get_dkz() { return m_header.dkz; };

    /**
    * \brief Get kx
    *
    * @param i Array Index in x-direction
    * @param t_i Translated array index in x-direction
    * @param k_x kx
    */
    void Get_kx( const ptrdiff_t i, ptrdiff_t& t_i, double& k_x )
    {
      t_i = (i+m_shift_x)%m_dimX;
      k_x = m_dkx*double(t_i-m_shift_x);
    }

    /**
    * \brief Get ky
    *
    * @param i Array Index in y-direction
    * @param t_i Translated array index in y-direction
    * @param k_y kyp
    */
    void Get_ky( const ptrdiff_t j, ptrdiff_t& t_j, double& k_y )
    {
      t_j = (j+m_loc_start_dimY+m_shift_y)%m_dimY;
      k_y = m_dky*double(t_j-m_shift_y);
    }

    /**
    * \brief Get kx
    *
    * @param i Array Index in x-direction
    */
    double Get_kx( const ptrdiff_t i )
    {
      return m_dkx*double(((i+m_shift_x)%m_dimX)-m_shift_x);
    }

    /**
    * \brief Get ky
    *
    * @param i Array Index in y-direction
    */
    double Get_ky( const ptrdiff_t j )
    {
      return m_dky*double(((j+m_loc_start_dimY+m_shift_y)%m_dimY)-m_shift_y);
    }

    /**
    * \brief Read Data from file
    *
    * @param filename Filename
    */
    void Read_File( std::string filename )
    {
      MPI_File   fh;
      ptrdiff_t  loc_offset;
      ptrdiff_t  loc_n;

      if( m_fs )
      {
        loc_offset = m_loc_start_dimY*m_dimX*m_dimZ;
        loc_n = m_loc_dimY*m_dimX*m_dimZ;
      }
      else
      {
        loc_offset = m_loc_start_dimX*m_dimY*m_dimZ;
        loc_n = m_loc_dimX*m_dimY*m_dimZ;
      }

      MPI_Offset offset = sizeof(generic_header) + sizeof(fftw_complex)*loc_offset;
      MPI_File_open( MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh );

      MPI_File_read_at( fh, offset, (double*)m_data, 2*loc_n, MPI_DOUBLE, MPI_STATUS_IGNORE );
      MPI_File_close( &fh );
    }

    /**
    * \brief Write Data to file
    *
    * @param filename Filename
    */
    void Write_File( std::string filename )
    {
      MPI_Status status;
      MPI_File   fh;
      ptrdiff_t  loc_offset;
      ptrdiff_t  loc_n;

      generic_header header=m_header;
      header.t = *m_time;

      if( m_fs )
      {
        loc_offset = m_loc_start_dimY*m_dimX*m_dimZ;
        loc_n = m_loc_dimY*m_dimX*m_dimZ;
        header.nDimX = m_header.nDimY;
        header.nDimY = m_header.nDimX;
        header.dx = m_header.dky;
        header.dy = m_header.dkx;
      }
      else
      {
        loc_offset = m_loc_start_dimX*m_dimY*m_dimZ;
        loc_n = m_loc_dimX*m_dimY*m_dimZ;
      }

      MPI_Offset offset = sizeof(generic_header) + sizeof(fftw_complex)*loc_offset;

      MPI_File_open( MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );

      if( m_rank == 0 )
        MPI_File_write( fh, &header, sizeof(generic_header), MPI_BYTE, &status );

      MPI_File_write_at( fh, offset, (double*)m_data, 2*loc_n, MPI_DOUBLE, MPI_STATUS_IGNORE );
      MPI_File_close( &fh );
    }
  protected:
    generic_header m_header; /// Header

    ptrdiff_t m_dimX; /// Global number of sampling points in x-direction
    ptrdiff_t m_dimY; /// Global number of sampling points in y-direction
    ptrdiff_t m_dimZ; /// Global number of sampling points in z-direction
    ptrdiff_t m_shift_x; /// Index shift in x-direction
    ptrdiff_t m_shift_y; /// Index shift in y-direction
    ptrdiff_t m_shift_z; /// Index shift in z-direction
    ptrdiff_t m_loc_dimX; /// Local number of sampling points in x-direction
    ptrdiff_t m_loc_dimY; /// Local number of sampling points in y-direction
    ptrdiff_t m_loc_start_dimX; /// Local Index Offset in x-direction
    ptrdiff_t m_loc_start_dimY; /// Local Index Offset in y-direction
    ptrdiff_t m_offset_rs; /// Offset in real space
    ptrdiff_t m_offset_fs; /// Offset in fourier space
    ptrdiff_t m_loc_n_rs; // Local number of elements in real space
    ptrdiff_t m_loc_n_fs; // Local number of elements in fourier space
    ptrdiff_t m_dimXZ; /// Product m_dimX * m_dimZ
    ptrdiff_t m_dimYZ; /// Product m_dimY * m_dimZ

    int m_rank; /// rank of the calling process in the group of comm
    int m_nprocs; /// number of processes in the group of comm

    bool m_fs; /// True if in fourier space

    double m_dx; /// Stepsize in x-direction
    double m_dy; /// Stepsize in y-direction
    double m_dz; /// Stepsize in z-direction
    double m_dkx; /// Stepsize in kx-direction
    double m_dky; /// Stepsize in ky-direction
    double m_dkz; /// Stepsize in kz-direction

    const double * m_time; /// Time

    fftw_complex * m_data; /// Data Array

    fftw_plan m_forwardPlan; /// FFTW Plan for forward fourier transformation
    fftw_plan m_backwardPlan; /// FFTW Plan for backward fourier transformation
  };
}}
