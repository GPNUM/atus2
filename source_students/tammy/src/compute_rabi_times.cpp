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
  table.clear();

  ifstream file( filename.c_str() );
  if( !file.is_open() )
  {
    cout << "Could not open file: " << filename << "." << endl;
    return false;
  }

  string line;
  getline(file, line);
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

  for( vector<Row>::iterator it = table.begin(); it != table.end();  )
  {
    Row row = *it;
    if( row.size() != 4 ) it = table.erase( it );
    else ++it;
  }

  return true;
}

void graph_walker( const unsigned dim, double * t, double * y, vector<double>& lo, vector<double>& hi, vector<double>& lo2, vector<double>& hi2 )
{
  const unsigned N = dim;
  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, dim);
  gsl_spline_init(spline, t, y, dim);

  double dt = (t[dim-1]-t[0])/double(N-1);

  cout << "dt == " << dt << endl;

  double abl_alt = gsl_spline_eval_deriv (spline, dt, acc);
  double abl_neu;

  for( unsigned i=2; i<N-2; i++ )
  {
    abl_neu = gsl_spline_eval_deriv (spline, double(i)*dt, acc);
    if( abl_neu*abl_alt < 0 ) 
    {
      if( abl_alt < 0 )
      {
        lo.push_back(double(i-5)*dt);  
        hi.push_back(double(i+4)*dt);
        i+=4;
      }
      else
      {
        lo2.push_back(double(i-5)*dt);  
        hi2.push_back(double(i+4)*dt);
        i+=4;
      }
    }
    abl_alt = abl_neu;
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);  
}

struct my_params
{
  gsl_interp_accel * acc;
  gsl_spline * spline;
};

double myfun_min ( double t, void * params )
{
  my_params * p = reinterpret_cast<my_params*>(params);
  return gsl_spline_eval (p->spline, t, p->acc);
}

double myfun_max ( double t, void * params )
{
  my_params * p = reinterpret_cast<my_params*>(params);
  return -gsl_spline_eval (p->spline, t, p->acc);
}

void find ( const unsigned dim, double * t, double * y, const vector<double>& lo, const vector<double>& hi, vector<double>& tau, bool mode ) 
{
  my_params parameter;

  int iter=0, max_iter=500, status;

  parameter.acc = gsl_interp_accel_alloc();
  parameter.spline = gsl_spline_alloc (gsl_interp_cspline, dim);
  gsl_spline_init(parameter.spline, t, y, dim);

  gsl_function F;
  if( mode )
    F.function = &myfun_min;
  else
    F.function = &myfun_max;
  F.params = &parameter;

  for( int i=0; i<lo.size(); i++ )
  {
    double a = lo[i];
    double b = hi[i];
    double m = 0.5*(lo[i]+hi[i]);

    // double fa = gsl_spline_eval (parameter.spline, a, parameter.acc);
    // double fm = gsl_spline_eval (parameter.spline, m, parameter.acc);
    // double fb = gsl_spline_eval (parameter.spline, b, parameter.acc);
    // printf( "f(%6.1g) = %10.5g | f(%6.1g) = %10.5g | f(%6.1g) = %10.5g\n", a, fa, m, fm, b, fb  ); 

    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, m, a, b);

    iter=0;
    do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      //printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, a, b, m, b - a);      

      status = gsl_min_test_interval (a, b, 1e-4, 0.0);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if (status == GSL_SUCCESS) tau.push_back(m);
    gsl_min_fminimizer_free (s);
  }

  gsl_spline_free (parameter.spline);
  gsl_interp_accel_free (parameter.acc);   
}

void partno ( const unsigned dim, double * t, double * y, const vector<double>& tau, vector<double>& f_tau ) 
{
  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, dim);
  gsl_spline_init(spline, t, y, dim);

  for( auto t : tau )
  {
    f_tau.push_back(gsl_spline_eval(spline, t, acc));
  }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc); 
}

int main(int argc, char *argv[])
{
  vector<Row> table;
  vector<double> lo, lo2, hi, hi2, tau, tau2, f_tau, f_tau2;

  bool bred = ReadFile( string(argv[1]), table );

  if( bred == false )
    return EXIT_FAILURE;
  
  const unsigned dim = table.size();

  cout << "dim == " << dim << endl;
  
  double *t = new double[dim];
  double *y = new double[dim];
  double *y2 = new double[dim];
 
  for( unsigned i=0; i<dim; i++ )
  {
    Row row = table[i];
    t[i] = row[0];
    y[i] = row[1];
    y2[i] = row[2];
  }

  // Impulszustand nullter Ordnung
  
  graph_walker( dim, t, y, lo, hi, lo2, hi2 );
  cout << "total number bracketing regions: " << hi.size() + hi2.size() << endl;  
  // for ( unsigned i=0; i<lo.size(); i++ )
  //   cout << "lo = " << lo[i] << ", hi = " << hi[i] << endl;    
  // for ( unsigned i=0; i<lo2.size(); i++ )
  //   cout << "lo2 = " << lo2[i] << ", hi2 = " << hi2[i] << endl;    

  tau.push_back(0);
  find( dim, t, y, lo, hi, tau, true );
  find( dim, t, y, lo2, hi2, tau, false );

  sort(tau.begin(), tau.end());
  cout << "total number of critical points " << tau.size() << endl;
  
  // cout << "critical points " << endl;
  // for ( unsigned i=0; i<tau.size(); i++ )
  //   cout << "tau = " << tau[i] << endl;    

  cout << "critical points diff" << endl;
  for ( unsigned i=1; i<tau.size(); i++ )
    cout << "tau = " << tau[i]-tau[i-1] << endl;

  partno (dim, t, y, tau, f_tau );         

  // Impulszustand erster Ordnung
  lo.clear();
  lo2.clear();
  hi.clear();
  hi2.clear();
  tau2.push_back(0);
  graph_walker( dim, t, y2, lo, hi, lo2, hi2 );
  find( dim, t, y2, lo, hi, tau2, true );
  find( dim, t, y2, lo2, hi2, tau2, false );

  sort(tau2.begin(), tau2.end());
  cout << "total number of critical points " << tau.size() << endl;
  
  // cout << "critical points " << endl;
  // for ( unsigned i=0; i<tau2.size(); i++ )
  //   cout << "tau2 = " << tau2[i] << endl;    

  cout << "critical points diff" << endl;
  for ( unsigned i=1; i<tau2.size(); i++ )
    cout << "tau2 = " << tau2[i]-tau2[i-1] << endl;    

  partno (dim, t, y, tau2, f_tau2 );         


  // Dateiausgabe
  string src_str(argv[1]);
  src_str.erase (src_str.begin()+9, src_str.end());
  
  ofstream outfile( src_str + "0.txt" );
  for( int i=1; i<tau.size(); i++ )
    outfile <<  tau[i]-tau[i-1] << endl;
  outfile.close();

  outfile.open( src_str + "1.txt" );
  for( int i=1; i<tau2.size(); i++ )
    outfile <<  tau2[i]-tau2[i-1] << endl;    
  outfile.close();

  outfile.open( src_str + "0_PartNo.txt" );
  for( int i=1; i<tau.size(); i++ )
    outfile <<  tau[i] << "\t" << f_tau[i] << endl;    
  outfile.close();

  outfile.open( src_str + "1_PartNo.txt" );
  for( int i=1; i<tau.size(); i++ )
    outfile <<  tau2[i] << "\t" << f_tau2[i] << endl;    
  outfile.close();

  delete [] t;
  delete [] y;
  delete [] y2;
return EXIT_SUCCESS;
}