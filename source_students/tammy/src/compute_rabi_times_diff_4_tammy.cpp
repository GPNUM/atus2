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

bool ReadFile( const std::string& filename, vector<vector<Row>>& data_tables )
{
  vector<Row> table;

  ifstream file( filename.c_str() );
  if( !file.is_open() )
  {
    cout << "Could not open file: " << filename << "." << endl;
    return false;
  }

  string line;
  vector<string> vec;

  getline(file, line); // skip the first line
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
  data_tables.push_back(table);
  return true;
}

int compute_differences( const vector<Row>& table_1, const vector<Row>& table_2, vector<vector<Row>>& dest, bool rel=false )
{
  // if( table_1.size() != table_2.size() ) return -1;

  vector<Row> new_table;
  Row new_row;

  if( rel )
  {
    for( int i=0; i<table_1.size(); i++ ) // rows
    {
      Row Row_table_1 = table_1[i];
      Row Row_table_2 = table_2[i];

      new_row.clear();
      new_row.push_back(Row_table_1[0]); // particle number of first species

      for( int j=2; j<=5; j++ ) // cols
      { 
        new_row.push_back( (Row_table_1[j]-Row_table_2[j])/Row_table_1[0] );
      }
      new_table.push_back(new_row);
    }
  }
  else
  {
    for( int i=0; i<table_1.size(); i++ ) // rows
    {
      Row Row_table_1 = table_1[i];
      Row Row_table_2 = table_2[i];

      new_row.clear();
      new_row.push_back(Row_table_1[0]); // particle number of first species

      for( int j=2; j<=5; j++ ) // cols
      { 
        new_row.push_back( Row_table_1[j]-Row_table_2[j] );
      }
      new_table.push_back(new_row);
    }
  }
  dest.push_back(new_table);
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
  }
}

int main(int argc, char *argv[])
{
  string filename;
  string basepath_1 = "Rb_K/";
  string basepath_2 = "Rb_K_nonlin_an/";
  vector<string> filenames;
  vector<vector<Row>> data;
  vector<vector<Row>> data_diff;

  const int ddetuning=500;

  for( int j=2; j<6; j++ )
  {
    filenames.clear();
    filenames.push_back(basepath_1 + "Rb_0__detuning_" + to_string(j*ddetuning) + "_tau.txt");
    filenames.push_back(basepath_2 + "Rb_0__detuning_" + to_string(j*ddetuning) + "_tau.txt");
    filenames.push_back(basepath_1 + "Rb_0__detuning_" + to_string(j*ddetuning) + "_N_of_tau.txt");
    filenames.push_back(basepath_2 + "Rb_0__detuning_" + to_string(j*ddetuning) + "_N_of_tau.txt");
    filenames.push_back(basepath_1 + "Rb_1__detuning_" + to_string(j*ddetuning) + "_tau.txt");
    filenames.push_back(basepath_2 + "Rb_1__detuning_" + to_string(j*ddetuning) + "_tau.txt");
    filenames.push_back(basepath_1 + "Rb_1__detuning_" + to_string(j*ddetuning) + "_N_of_tau.txt");
    filenames.push_back(basepath_2 + "Rb_1__detuning_" + to_string(j*ddetuning) + "_N_of_tau.txt");

    for( auto name : filenames )
    {
      if(!ReadFile( name, data ))
      {
        cout << "unable to open file " << name << endl;
        return EXIT_SUCCESS;
      }
    }

    compute_differences( data[0], data[1], data_diff );   
    compute_differences( data[2], data[3], data_diff, true );   
    compute_differences( data[4], data[5], data_diff );   
    compute_differences( data[6], data[7], data_diff, true );   

    output_table( "Rb_0__detuning_" + to_string(j*ddetuning) + "_tau_diff.txt", data_diff[0] );   
    output_table( "Rb_0__detuning_" + to_string(j*ddetuning) + "_N_of_tau_diff.txt", data_diff[1] );   
    output_table( "Rb_1__detuning_" + to_string(j*ddetuning) + "_tau_diff.txt", data_diff[2] );   
    output_table( "Rb_1__detuning_" + to_string(j*ddetuning) + "_N_of_tau_diff.txt", data_diff[3] );   
  }
return EXIT_SUCCESS;
}