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


/** @file
  *
  * \class CRT_shared_mpi
  * Class with basic definitions and functions in MPI.
  * Supports only two and three dimensions.
  * Global size is the total number of points.
  * Local size is the number of points of the process.
  */
#ifndef __class_CRT_shared_mpi__
#define __class_CRT_shared_mpi__

#include "my_structs.h"
#include "CPoint.h"
#include "fftw3-mpi.h"
#include <cmath>
#include <fstream>
#include <cassert>

/** Function pointer with a sequence_item
  */
typedef void (*StepFunction)(void *,sequence_item &);

class CRT_shared_mpi
{
public:
  /** Constructor. Gives the rank of the process in communicator group MPI_COMM_WORLD */
  CRT_shared_mpi() :
    m_header( {}),
            m_dimX(0),
            m_dimY(0),
            m_dimZ(0),
            m_loc_dimX(0),
            m_loc_dimY(0),
            m_loc_start_dimX(0),
            m_loc_start_dimY(0),
            m_no_of_pts(0),
            m_no_of_pts_fs(0),
            m_shift_x(0),
            m_shift_y(0),
            m_shift_z(0),
            m_alloc(0),
            m_ar(0),
            m_ar_k(0)
  {
    MPI_Comm_rank( MPI_COMM_WORLD, &m_myrank );
  };

  double &Get_t()
  {
    return m_header.t;
  };
  const double &Get_dt() const
  {
    return m_header.dt;
  };
  const double &Get_dx() const
  {
    return m_header.dx;
  };
  const double &Get_dy() const
  {
    return m_header.dy;
  };
  const double &Get_dz() const
  {
    return m_header.dz;
  };
  const double &Get_dkx() const
  {
    return m_header.dkx;
  };
  const double &Get_dky() const
  {
    return m_header.dky;
  };
  const double &Get_dkz() const
  {
    return m_header.dkz;
  };

  /// Get number of grid points in x-direction
  long long &Get_dimX()
  {
    return m_header.nDimX;
  };
  long long &Get_dimY()
  {
    return m_header.nDimY;
  };
  long long &Get_dimZ()
  {
    return m_header.nDimZ;
  };
  /// Get number of local points
  ptrdiff_t Get_No_Points()
  {
    return m_no_of_pts;
  };

  /// Change size of time steps to dt and call Init() to recalculate the kinetic operator
  void Set_dt( const double dt )
  {
    m_header.dt = dt;
    Init();
  };
protected:
  /** Read header of a file into #m_header
    *
    * @param filename Name of the file to be read
    * @param dim Dimensions of the file
    */
  void Read_header(const std::string &filename, const int dim)
  {
    //Read header into m_header
    std::ifstream file1( filename, std::ifstream::binary );
    file1.read( (char *)&m_header, sizeof(generic_header) );
    file1.close();

    switch ( dim )
    {
    //2D
    case 2:
      m_ar = m_header.dx*m_header.dy;
      m_ar_k = m_header.dkx*m_header.dky;
      // Get local data size
      m_alloc = fftw_mpi_local_size_2d_transposed( Get_dimX(), Get_dimY(), MPI_COMM_WORLD, &m_loc_dimX, &m_loc_start_dimX, &m_loc_dimY, &m_loc_start_dimY );
      m_no_of_pts = m_loc_dimX*m_header.nDimY;
      m_no_of_pts_fs = m_loc_dimY*m_header.nDimX;
      m_shift_x = m_header.nDimX/2;
      m_shift_y = m_header.nDimY/2;
      m_header.nDimZ = 1;
      break;
    //3D
    case 3:
      m_ar = m_header.dx*m_header.dy*m_header.dz;
      m_ar_k = m_header.dkx*m_header.dky*m_header.dkz;
      // Get local data size
      m_alloc = fftw_mpi_local_size_3d_transposed( Get_dimX(), Get_dimY(), Get_dimZ(), MPI_COMM_WORLD, &m_loc_dimX, &m_loc_start_dimX, &m_loc_dimY, &m_loc_start_dimY );
      m_no_of_pts = m_loc_dimX*Get_dimY()*Get_dimZ();
      m_no_of_pts_fs = m_loc_dimY*Get_dimX()*Get_dimZ();
      m_shift_x = m_header.nDimX/2;
      m_shift_y = m_header.nDimY/2;
      m_shift_z = m_header.nDimZ/2;
      break;
    }
  }

  /// The header of a file is read into this struct
  generic_header m_header;
  /// Global number of points in x-direction (deprecated?)
  ptrdiff_t m_dimX;
  ptrdiff_t m_dimY;
  ptrdiff_t m_dimZ;
  /// Local number of points in x-direction
  ptrdiff_t m_loc_dimX;
  ptrdiff_t m_loc_dimY;
  /// Local grid starts at position m_loc_start_dimX of global grid
  ptrdiff_t m_loc_start_dimX;
  ptrdiff_t m_loc_start_dimY;
  /// Number of local points
  ptrdiff_t m_no_of_pts;
  /// number of local points in fourier space
  ptrdiff_t m_no_of_pts_fs;
  /// Shift of x-axis
  ptrdiff_t m_shift_x;
  ptrdiff_t m_shift_y;
  ptrdiff_t m_shift_z;
  // Local data size
  ptrdiff_t m_alloc;

  /// Rank of MPI process
  int m_myrank;

  /// Volume element
  double m_ar;
  /// Volume element in Fourierspace
  double m_ar_k;
  ///Calculates the kinetic operator
  virtual void Init()=0;
};
#endif
