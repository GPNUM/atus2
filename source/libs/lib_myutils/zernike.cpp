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


#include "zernike.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <boost/math/special_functions/binomial.hpp>

using namespace boost::math;

zernike::zernike( const generic_header &header )
{
}

void zernike::rescale_zern( double scale)
{
  for ( int i = 0; i<m_zern.size(); i++)
  {
    for ( int j = 0; j<m_zern[i].size(); j++)
    {
      m_zern[i][j] *= scale;
    }
  }
}

void zernike::set_max( double val)
{

}

void zernike::add_zern( int n, int m, double scaling )
{
  int p, q, l, xpow, ypow, xy;
  l = n-2*m;

  double factor;

  if ( l <= 0 )
  {
    p = 0;
    q = ( n%2 == 0 ? -l/2 : (-l-1)/2 );
  }
  else
  {
    p = 1;
    q = ( n%2 == 0 ? l/2 - 1: (l-1)/2 );
  }

  l = (l<0 ? -l : l);
  m =(n-l)/2;

  for ( int i=0; i<=q; i++ )
  {
    for ( int j=0; j<=m; j++ )
    {
      for ( int k=0; k<=(m-j); k++ )
      {
        factor = ( (i+j)%2 == 0 ? 1 : -1 );
        factor *= binomial_coefficient<double>( l, 2*i+p );
        factor *= binomial_coefficient<double>( m-j, k );
        factor *= factorial<double>(n-j)/(factorial<double>(m-j)*\
                                          factorial<double>(m-j)*factorial<double>(n-m-j));
        ypow = 2 * (i+k) + p;
        xpow = n - 2 * (i+j+k) - p;

        resize_zern( xpow+1, ypow +1);
        m_zern[xpow][ypow] += factor*scaling;
      }
    }
  }
}

void zernike::read_zernike()
{
  ifstream f("zernike.txt");
  string line;

  while (getline(f, line))
  {
    istringstream ss(line);

    int n, m;
    double scaling;

    ss >> n >> m >> scaling;
    add_zern(n,m,scaling);
  }

}

void zernike::resize_zern( int x, int y )
{
  int element;
  x >= y ? element = x : element = y;
  if ( element > m_zern.size() )
  {
    m_zern.resize(element);
    for ( int i = 0; i< element; ++i)
      m_zern[i].resize(element);
  }
}

//void zernike::Get_zern_coord( int x, int y int$ element)
//{
//  int m = (x+y) + 1;
//  m = m*(m+1)/2; //Position (x,0) Element
//  element = m + y; //Position (x,y) Element

//  if( element > m_zern.size() )
//    m_zern.resize(element);
//}

//void zernike::Get_power( int& x, int& y, int element)
//{

//}

void zernike::generate_zern()
{
}

void zernike::calc_zernike( const generic_header &header, double *field )
{
  int N = header.nDimX;

  for ( int i = 0; i < N; i++)
  {
    for ( int j = 0; j < N; j++)
    {
      int ij = j+N*i;
      double x = (i-N/2)*header.dx/header.xMax;
      double y = (j-N/2)*header.dy/header.yMax;
      for ( int xpow = 0; xpow<m_zern.size(); xpow++)
      {
        for ( int ypow = 0; ypow<m_zern[xpow].size(); ypow++)
        {
          field[ij] += m_zern[xpow][ypow]*(pow(x,xpow)*pow(y,ypow));
        }
      }
    }
  }
}



void zernike::print_zernike()
{
  for ( int i = 0; i<m_zern.size(); i++)
  {
    for ( int j = 0; j<m_zern[i].size(); j++)
    {
      std::cout << m_zern[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
