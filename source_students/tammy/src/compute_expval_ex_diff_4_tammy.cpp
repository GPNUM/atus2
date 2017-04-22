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

/*
Želimir Marojević
*/

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "strtk.hpp"
 
using namespace std; 

typedef vector<double> Row;

bool ReadFile( const std::string& filename, vector<Row>& table )
{
  ifstream file( filename.c_str() );
  if( !file.is_open() )
  {
    cout << "Could not open file: " << filename << "." << endl;
    return false;
  }

  string line;
  vector<string> vec;

  while(file)
  {
    line.clear();
    vec.clear();
    getline(file, line);
    strtk::parse(line,"\t",vec);
    if( line.size() == 0 ) continue;

    Row row;
    double data;
    for( auto i : vec )
    {
      try
      { 
        data = stod(i);
      }
      catch( const std::invalid_argument& ia )
      {
        continue;
      }      
      row.push_back(data);
    }
    table.push_back(row);
  }
  return true;
}

int compute_differences( const vector<Row>& table_1, const vector<Row>& table_2, vector<Row>& new_table )
{
  // if( table_1.size() != table_2.size() ) return -1;
  Row new_row;

  for( int i=0; i<table_1.size(); i++ ) // rows
  {
    Row Row_table_1 = table_1[i];
    Row Row_table_2 = table_2[i];

    new_row.clear();
    new_row.push_back(Row_table_1[0]); // particle number of first species
    new_row.push_back(Row_table_1[1]); // detuning

    for( int j=2; j<Row_table_1.size(); j++ ) // cols
    { 
      new_row.push_back( Row_table_1[j]-Row_table_2[j] );
    }
    new_table.push_back(new_row);
  }
return 0;
}

void output_table( const string& filename, const vector<Row>& data )
{
  ofstream outfile(filename);

  for( int i=0; i<data.size(); i++ )
  {
    Row row_i = data[i];
    outfile << row_i[0];
    for( int j=1; j<row_i.size(); j++ )
    {
      outfile << "\t" << row_i[j];
    }
    outfile << endl;
    //if( (i+1) % 17 == 0 ) outfile << endl;
  }
}

int main(int argc, char *argv[])
{
  if( argc != 3 )
  {
    cout << "no filenames provided" << endl;
    return EXIT_SUCCESS;
  }

  string filename_1 = argv[1];
  string filename_2 = argv[2];

  vector<Row> table_1, table_2, table_diff;
  vector<Row> data_;

  if(!ReadFile( filename_1, table_1 ))
  {
    cout << "unable to open file " << filename_1 << endl;
    return EXIT_SUCCESS;
  }

  if(!ReadFile( filename_2, table_2 ))
  {
    cout << "unable to open file " << filename_2 << endl;
    return EXIT_SUCCESS;
  }

  compute_differences( table_1, table_2, table_diff );   

  output_table( "difference.txt", table_diff );   
return EXIT_SUCCESS;
}