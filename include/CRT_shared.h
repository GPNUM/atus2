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


/** @file */

#ifndef __class_CRT_shared__
#define __class_CRT_shared__

#include "my_structs.h"
#include "CPoint.h"
#include "fftw3.h"
#include <cmath>
#include <fstream>
#include <cassert>
#include <array>

/** Function pointer with a sequence_item
  */
typedef void (*StepFunction)(void *,sequence_item &);

/** Class with basic definitions and functions for one, two or three dimensions
  *
  */
class CRT_shared
{
public:
  /** Constructor
    * Open log file
    */
  CRT_shared() :
    m_no_of_pts(0),
    m_no_of_pts_red(0),
    m_shift_x(0),
    m_shift_y(0),
    m_shift_z(0),
    m_ar(0),
    m_ar_k(0)
  {
    m_header = {};
    m_log.open("log.txt");
    m_log << std::scientific;
    m_log << std::setprecision(10);
  }

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
  /// Get number of grid points in x direction
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
  /// Get total number of grid points
  const int &Get_No_Points() const
  {
    return m_no_of_pts;
  };

  /// Change size of time steps to dt and call Init() to recalculate the kinetic operator
  void Set_dt( const double dt )
  {
    m_header.dt = dt;
    Init();
  };
  ///log file
  std::ofstream m_log;

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
    if ( !file1.is_open() ) throw std::string( "Could not open file " + filename + ".\n" );
    file1.read( (char *)&m_header, sizeof(generic_header) );
    file1.close();

    switch ( dim )
    {
    //1D
    case 1:
      m_ar = m_header.dx;
      m_ar_k = m_header.dkx;
      m_no_of_pts = m_header.nDimX;
      m_shift_x = m_header.nDimX/2;
      m_no_of_pts_red = m_shift_x+1;
      m_header.nDimY = 1;
      m_header.nDimZ = 1;
      break;
    //2D
    case 2:
      m_ar = m_header.dx*m_header.dy;
      m_ar_k = m_header.dkx*m_header.dky;
      m_no_of_pts = m_header.nDimX*m_header.nDimY;
      m_shift_x = m_header.nDimX/2;
      m_shift_y = m_header.nDimY/2;
      m_no_of_pts_red = m_header.nDimX*(m_shift_y+1);
      m_header.nDimZ = 1;
      break;
    //3D
    case 3:
      m_ar = m_header.dx*m_header.dy*m_header.dz;
      m_ar_k = m_header.dkx*m_header.dky*m_header.dkz;
      m_no_of_pts = m_header.nDimX*m_header.nDimY*m_header.nDimZ;
      m_shift_x = m_header.nDimX/2;
      m_shift_y = m_header.nDimY/2;
      m_shift_z = m_header.nDimZ/2;
      m_no_of_pts_red = m_header.nDimX*m_header.nDimY*(m_shift_z+1);
      break;
    }
  }

  /// The header of a file is read into this struct
  generic_header m_header;
  /// Total number of points
  int m_no_of_pts;
  int m_no_of_pts_red;
  /// Shift of x-axis
  int m_shift_x;
  int m_shift_y;
  int m_shift_z;
  /// Volume element
  double m_ar;
  /// Volume element in Fourierspace
  double m_ar_k;
  ///Calculates kinetic operator
  virtual void Init()=0;
};
#endif
