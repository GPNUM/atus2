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

#include <stdio.h>
#include <iostream>
#include "cmath"
#include "lis.h"
#include <fstream>
#include "fftw3.h"
#include "my_structs.h"

using namespace std;

class lis_sg
{
  public:
    lis_sg();
    ~lis_sg();

    double Particle_number( LIS_VECTOR& );
    double coeff_a( const double x, const double t );

    void setup_initial_wavefunction( LIS_VECTOR& );
    void assemble_system_and_rhs( LIS_VECTOR& );
    void set_boundary_values( LIS_VECTOR& );
    void solve( LIS_VECTOR& );
    void output_vector( const std::string&, LIS_VECTOR& );
    void run( const double );

  protected:
    unsigned m_dim;

    double m_t;
    double m_dt;
    double m_dx;
    double m_xmin;
    double m_xmax;
    double m_alpha;

    LIS_MATRIX m_Dop;
    LIS_MATRIX m_Dop_2;
    LIS_VECTOR m_Psi_1;
    LIS_VECTOR m_Psi_2;
    LIS_VECTOR m_rhs;
    LIS_SOLVER m_solver;
};

lis_sg::lis_sg()
{
  m_dim = 1024;
  m_xmin = -10;
  m_xmax = 10;

  m_dx = (m_xmax - m_xmin)/double(m_dim-1);
  m_dt = 0.01;
  m_t = 0;

  m_alpha = 1;

  lis_vector_create(0, &m_Psi_1);
  lis_vector_set_size(m_Psi_1, 0, 2*m_dim);

  lis_vector_create(0, &m_Psi_2);
  lis_vector_set_size(m_Psi_2, 0, 2*m_dim);

  lis_vector_create(0, &m_rhs);
  lis_vector_set_size(m_rhs, 0, 2*m_dim);

  lis_solver_create(&m_solver);
  //lis_solver_set_option("-print all",m_solver);
  lis_solver_set_option("-i GMRES",m_solver);
  //lis_solver_set_option("-p ilu",m_solver);
  lis_solver_set_option("-maxiter 10000",m_solver);

  setup_initial_wavefunction(m_Psi_1);
}

lis_sg::~lis_sg()
{
  lis_matrix_destroy(m_Dop);
  lis_vector_destroy(m_Psi_1);
  lis_vector_destroy(m_Psi_2);
  lis_vector_destroy(m_rhs);
  lis_solver_destroy(m_solver);
  lis_finalize();
}

double lis_sg::Particle_number( LIS_VECTOR& vec )
{
  double retval=0, re, im;

  for( int i=0; i<m_dim; i+=2 )
  {
    lis_vector_get_value(vec,2*i,&re);
    lis_vector_get_value(vec,2*i+1,&im);
    retval += re*re + im*im;
  }
return retval*m_dx;
}

void lis_sg::set_boundary_values( LIS_VECTOR& b )
{
  lis_vector_set_value(LIS_INS_VALUE, 0, 0.0, b);
  lis_vector_set_value(LIS_INS_VALUE, 1, 0.0, b);
  lis_vector_set_value(LIS_INS_VALUE, 2*m_dim-2, 0.0, b);
  lis_vector_set_value(LIS_INS_VALUE, 2*m_dim-1, 0.0, b);
}

void lis_sg::setup_initial_wavefunction( LIS_VECTOR& x )
{
  for( int i=0; i<m_dim; i++ )
  {
    double xval = m_xmin + double(i)*m_dx;
    double re = 1.0/sqrt(2*3.1415926535897932)*exp(-(xval*xval)/2.0);
    double im = 0.0;
    lis_vector_set_value(LIS_INS_VALUE, 2*i, re, x);
    lis_vector_set_value(LIS_INS_VALUE, 2*i+1, im, x);
  }
}

double lis_sg::coeff_a( const double x, const double t )
{
return sin(0.1*t)*sin(x-0.01*t);
}

void lis_sg::assemble_system_and_rhs( LIS_VECTOR& rhs )
{
  int N = m_dim/2;
/*
  const double diag = 1;
  const double fak0 = 30*m_alpha*m_dt/(12*m_dx*m_dx);
  const double fak1 = 16*m_alpha*m_dt/(12*m_dx*m_dx);
  const double fak2 = m_alpha*m_dt/(12*m_dx*m_dx);
*/
  const double diag = m_dx*m_dx/m_dt;
  const double fak0 = 30.0/12.0;
  const double fak1 = 16.0/12.0;
  const double fak2 = 1.0/12.0;

  double x = m_xmin;
  double a = coeff_a( x, m_t );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0, diag, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0+1, -a*fak0, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0+3, a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0+5, -a*fak2, m_Dop);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1-1, -a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1, diag, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1+1, a*fak0, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1+3, -a*fak2, m_Dop);

  x += m_dx;
  a = coeff_a( x, m_t );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2-1, a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2  , diag, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2+1, -a*fak0, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2+3, a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2+5, -a*fak2, m_Dop);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)-3, -a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)-1, a*fak0, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)  , diag, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)+1, -a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)+3, a*fak2, m_Dop);

  for ( int i = 2; i < (m_dim-2); i++ ) 
  {
    x += m_dx;
    a = coeff_a( x, m_t );
    // real
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)-3, -a*fak2, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)-1, a*fak1, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)  , diag, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+1, -a*fak0, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+3, a*fak1, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+5, -a*fak2, m_Dop);
    // imag
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-5, a*fak2, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-3, -a*fak1, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-1, a*fak0, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)  , diag, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)+1, -a*fak1, m_Dop);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)+3, a*fak2, m_Dop);
  }

  x += m_dx;
  a = coeff_a( x, m_t+m_dt );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))-3, -a*fak2, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))-1, a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))  , diag, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))+1, -a*fak0, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))+3, a*fak1, m_Dop);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-5, a*fak2, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-3, -a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-1, a*fak0, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)  , diag, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)+1, -a*fak1, m_Dop);

  x += m_dx;
  a = coeff_a( x, m_t+m_dt );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))-3, -a*fak2, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))-1, a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1)), diag, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))+1, -a*fak0, m_Dop);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-5, a*fak2, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-3, -a*fak1, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-1, a*fak0, m_Dop);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1), diag, m_Dop);
  
  lis_matrix_set_type(m_Dop, LIS_MATRIX_CSR);
  lis_matrix_assemble(m_Dop);

  x = m_xmin;
  a = coeff_a( x, m_t+m_dt );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0, diag, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0+1, -a*fak0, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0+3, a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 0, 0+5, -a*fak2, m_Dop_2);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1-1, -a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1, diag, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1+1, a*fak0, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 1, 1+3, -a*fak2, m_Dop_2);

  x += m_dx;
  a = coeff_a( x, m_t+m_dt );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2-1, a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2  , diag, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2+1, -a*fak0, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2+3, a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2, 2+5, -a*fak2, m_Dop_2);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)-3, -a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)-1, a*fak0, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)  , diag, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)+1, -a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 3, (3)+3, a*fak2, m_Dop_2);

  for ( int i = 2; i < (m_dim-2); i++ ) 
  {
    x += m_dx;
    a = coeff_a( x, m_t+m_dt );
    // real
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)-3, -a*fak2, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)-1, a*fak1, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)  , diag, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+1, -a*fak0, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+3, a*fak1, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i, (2*i)+5, -a*fak2, m_Dop_2);
    // imag
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-5, a*fak2, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-3, -a*fak1, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)-1, a*fak0, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)  , diag, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)+1, -a*fak1, m_Dop_2);
    lis_matrix_set_value(LIS_INS_VALUE, 2*i+1, (2*i+1)+3, a*fak2, m_Dop_2);
  }

  x += m_dx;
  a = coeff_a( x, m_t+m_dt );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))-3, -a*fak2, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))-1, a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))  , diag, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))+1, -a*fak0, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2), (2*(N-2))+3, a*fak1, m_Dop_2);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-5, a*fak2, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-3, -a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)-1, a*fak0, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)  , diag, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, 2*(N-2)+1, (2*(N-2)+1)+1, -a*fak1, m_Dop_2);

  x += m_dx;
  a = coeff_a( x, m_t+m_dt );
  // real
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))-3, -a*fak2, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))-1, a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1)), diag, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)), (2*(N-1))+1, -a*fak0, m_Dop_2);
  // imag
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-5, a*fak2, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-3, -a*fak1, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1)-1, a*fak0, m_Dop_2);
  lis_matrix_set_value(LIS_INS_VALUE, (2*(N-1)+1), (2*(N-1)+1), diag, m_Dop_2);
  
  lis_matrix_set_type(m_Dop_2, LIS_MATRIX_CSR);
  lis_matrix_assemble(m_Dop_2);


  lis_matvect( m_Dop_2, rhs, m_rhs);
  set_boundary_values( m_rhs );

  //string filename = "matrix-" + to_string(m_t) + ".txt";
  //lis_output_matrix(m_Dop, LIS_FMT_MM, const_cast<char *>(filename.c_str()));
}

void lis_sg::solve( LIS_VECTOR& sol )
{

  LIS_INT err=lis_solve(m_Dop, m_rhs, sol, m_solver); CHKERR(err);

  // lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time);
  // lis_solver_get_residualnorm(solver,&resid);  //lis_solver_get_iterex(solver,&iter,&iter_double,&iter_quad);
  //  if (iter > 1000) {
  //    printf("Lis iter %i\n", iter);
}

void lis_sg::output_vector( const std::string& filename, LIS_VECTOR& vec )
{
  std::ofstream ofs(filename);
  for( int i=0; i<m_dim; i++ )
  {
    double re, im;
    lis_vector_get_value(vec,2*i,&re);
    lis_vector_get_value(vec,2*i+1,&im);
    ofs << m_xmin+double(i)*m_dx << "\t" << re*re+im*im << "\t" << re << "\t" << im << std::endl;
  }
}

void lis_sg::run(const double T)
{
   unsigned NT = unsigned(T/m_dx/0.5);

   double N = Particle_number( m_Psi_1 );
   cout << "N == " << N << endl;
 
   for( unsigned c=0; c < NT; c+=2 )
   {
     lis_matrix_create(0,&m_Dop);
     lis_matrix_set_size(m_Dop,0,2*m_dim);
     lis_matrix_create(0,&m_Dop_2);
     lis_matrix_set_size(m_Dop_2,0,2*m_dim);
     assemble_system_and_rhs( m_Psi_1 );
     solve( m_Psi_2 );
     lis_matrix_destroy(m_Dop);
     lis_matrix_destroy(m_Dop_2);

     lis_matrix_create(0,&m_Dop);
     lis_matrix_set_size(m_Dop,0,2*m_dim);
     lis_matrix_create(0,&m_Dop_2);
     lis_matrix_set_size(m_Dop_2,0,2*m_dim);
     assemble_system_and_rhs( m_Psi_2 );
     solve( m_Psi_1 );
     lis_matrix_destroy(m_Dop);
     lis_matrix_destroy(m_Dop_2);

     double N = Particle_number( m_Psi_1 );
     cout << "N == " << N << endl;
   }
}

int main(int argc, char *argv[])
{
  //LIS_INT err,iter,iter_double,iter_quad;
  //double time,itime,ptime,p_c_time,p_i_time;
  //LIS_REAL resid;

  lis_initialize(&argc, &argv);

  lis_sg sol;
  sol.run(10);

  return 0;
}
