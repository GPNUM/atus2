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

#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

int main ( int argc, char *argv[] )
{
  if ( argc != 4 )
  {
    cout << argv[0] <<  " ti dt tf" << endl;
    return EXIT_FAILURE;
  }

  vector<string> args;
  for ( int i=1; i<argc; i++ )
    args.push_back(argv[i]);

  const int ti = stod(args[0]);
  const int dt = stod(args[1]);
  const int tf = stod(args[2]);

  ofstream out("solution.pvd");
  out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n<Collection>\n";
  for ( int t=ti; t<=tf; t+=dt )
    out << "<DataSet timestep=\"" << to_string(t) << ".0" << "\" group=\"\" part=\"0\" file=\"3d_" << to_string(t) << ".000_1.vtk.vti\"/>\n";

  out << "</Collection>\n</VTKFile>\n";
  return EXIT_SUCCESS;
}