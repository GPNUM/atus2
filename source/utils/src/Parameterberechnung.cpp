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
#include <cmath>
#include <complex.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <strings.h>
#include <sys/time.h>
#include <omp.h>
#include "fftw3.h"
#include "my_structs.h"
#include "local_params.h"
#include "my_math_consts.h"
//#include "tinyxml.h"

using namespace std;

extern void Set_filename_2( char[], char[], TiXmlNode * );
extern void Set_localparams( local_params *, TiXmlNode * );
extern void Set_Stepsize( double &, TiXmlNode * );
extern void Set_EPS( double &, TiXmlNode * );


const int MAXITER  = 1000;
const double error = 1E-20;


int NewtonMethod_Dipol(double &x , double eta, double delta, double beta)
{
  int n = 1;
#define F(x) ( eta*delta*2*x*exp(-delta*x*x) + beta )
#define dF(x) ( eta*delta*2*exp(-delta*x*x) - eta*delta*delta*4*x*x*exp(-delta*x*x) )
  printf("Start: %E\n", x);

  while ( ( fabs(F(x)) > error ) && ( n <= MAXITER ) )
  {
    x = x - ( F(x) / dF(x) );
    printf("\t%E\n", x);
    n++;
  }
  printf("The end solution is: %E\n", x);
  return n;
}

int NewtonMethod_ham_osci(double &x , double om_xq, double beta)
{
  int n = 1;
#define G(x) ( om_xq*2*x + beta )
#define dG(x) ( om_xq*2 )
  printf("Start: %E\n", x);

  while ( ( fabs(G(x)) > error ) && ( n <= MAXITER ) )
  {
    x = x - ( G(x) / dG(x) );
    printf("\t%E\n", x);
    n++;
  }
  printf("The end solution is: %E\n", x);
  return n;
}

//-----------------------------------------------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{
  if ( argc != 2 )
  {
    printf( "No parameter xml file specified.\n\nGenerating example file: INPUT-EXAMPLE.xlm\n\n" );

    //Generate Input.xml-file
    ofstream outf_in( "INPUT_EXAMPLE.xml" );

    outf_in << setprecision(10);
    outf_in << "<SET_PARAMS>\n\n";
    outf_in << "<DATEI_NAME>3d_0.000_1.bin</DATEI_NAME>\n";
    outf_in << "<DATEI_NAME_2>3d_0.000_2.bin</DATEI_NAME_2>\n\n";
    outf_in << "<!-- First species: val=\"1\", second species: val=\"2\" or both species: val=\"12\" ? -->\n";
    outf_in << "<SPECIES val=\"12\"></SPECIES>\n\n";
    outf_in << "<!-- Gravitational acceleration (in units of g) -->\n";
    outf_in << "<GRAV val=\"1\"></GRAV>\n\n";
    outf_in << "<!-- LASERSYSTEM -->\n";
    outf_in << "<W0 val=\"0.00005\" ></W0>\n";
    outf_in << "<LAMBDA val=\"0.0000019\" ></LAMBDA>\n";
    outf_in << "<P val=\"0.3\" ></P>\n\n";
    outf_in << "<!-- Number of particles (_1: Rb, _2, K)  -->\n";
    outf_in << "<N val=\"330\" val_2=\"330\"></N>\n\n";
    outf_in << "<!-- BECs (_11: Rb-Rb, _12: Rb-K, _22: K-K)  -->\n";
    outf_in << "<A_A0 val=\"90.0\" val_12=\"28.0\" val_21=\"28.0\" val_2=\"-33.0\"></A_A0>\n\n";
    outf_in << "<!-- Atomic mass: val=m_Rb , val_2=m_K -->\n";
    outf_in << "<!-- in kg: m_Rb=1.443161930E-25 kg    m_K=6.470080165E-26 kg -->\n";
    outf_in << "<!-- in atomic mass unit (u=1.660538921e-27kg): -->\n";
    outf_in << "<M val=\"86.90925046977566\" val_2=\"38.96373691201183\"></M>\n\n";
    outf_in << "<!-- Resonant wavelength of the first exited state:  -->\n";
    outf_in << "<R val=\"0.000000780\" val_2=\"0.0000007665\"></R>\n\n";
    outf_in << "<!-- TRAP PARAMETERS (in case of harmonic trap) -->\n";
    outf_in << "<OM val_x=\"1111.0\" val_y=\"1111.0\" val_z=\"43.9822971\"></OM>\n\n";
    outf_in << "<!-- SIMULATION -->\n";
    outf_in << "<NK val=\"100\"></NK><!-- After NK steps an output file is written -->\n";
    outf_in << "<NA val=\"100\"></NA><!-- Output interval  -->\n";
    outf_in << "<MXITER val=\"400000\"></MXITER><!-- No of iterations -->\n";
    outf_in << "<EPS val=\"0.0001\"></EPS><!-- L2-Gradient error -->\n";
    outf_in << "<STEP_SIZE val=\"0.0001\"></STEP_SIZE><!-- Initial stepsize -->\n";
    outf_in << "<FIN_N val=\"1\" val_2=\"1\"></FIN_N>\n\n";
    outf_in << "<!-- Typical evolution time of BEC (depends of trap frequencies) -->\n";
    outf_in << "<T val=\"0.001\" ></T>\n";
    outf_in << "<!-- Start/end time of interaction with external potential (in time units given by T)  -->\n";
    outf_in << "<T_ON val=\"0.0\"></T_ON>\n";
    outf_in << "<T_OFF val=\"10.0\"></T_OFF>\n";
    outf_in << "<!-- Typical length scale -->\n";
    outf_in << "<L val=\"0.000001\"></L>\n\n";
    outf_in << "</SET_PARAMS>" << "\n";
    return 0;
  }

  local_params   params;
  char filename[100], filename2[100];
  int  lc = 0;
  double stepsize, eps;

  const double pi   = 3.14159265358979323846264338328;
  const double c    = 299792458;
  const double hbar = 1.054571596E-34;
  const double a_0  = 5.291772081e-11; //Bohrradius
  const double epsilon0 = 8.854187815E-12;
  const double ee = 1.602176462E-19;
  const double me = 9.109381882E-31;
  const double u = 1.660538921e-27; //atomic mass unit

  TiXmlDocument doc( argv[1] );
  doc.LoadFile();

  TiXmlNode *pChild = 0;
  for ( pChild = doc.FirstChild( "SET_PARAMS" ); pChild != 0; pChild = pChild->NextSibling() )
  {
    lc++;

    //User input: Parameter setting of simulation
    bzero( &params, sizeof(local_params) );
    Set_filename_2( filename, filename2, pChild );
    Set_localparams( &params, pChild );
    Set_Stepsize( stepsize, pChild );
    Set_EPS( eps, pChild );

    int Na = params.Na;
    int Nk = params.Nk;
    int MaxIter = params.MaxIter;
    int grav = params.grav;
    double w_0 = params.w_0;
    double P = params.P;
    double lambda = params.lambda;
    double o_m_x = params.o_m_x;
    double o_m_y = params.o_m_y;
    double o_m_z = params.o_m_z;
    double T = params.T;
    double t_off = params.t_off;
    double t_on = params.t_on;
    double L = params.l;
    double a_11 = params.a_a0*a_0;
    double a_22 = params.a_a0_2*a_0;
    double a_12 = params.a_a0_12*a_0;
    double a_21 = params.a_a0_21*a_0;
    int N_1 = params.Particle_Number;
    int N_2 = params.Particle_Number_2;
    double g_F  = 9.80665*grav;
    double m_1 = params.m_1*u;  // [kg]
    double m_2 = params.m_2*u;  // [kg]
    double r_1 = params.r_1;
    double r_2 = params.r_2;
    double species = params.species;

    //Laser system, i.e. experimental setup
    double omega = 2*pi*c/params.lambda;
    double z_R = (pi*w_0*w_0)/lambda ;
    double x, y, z;


//---------------------------------------------------------------------------
//Calculate nonlinear interaction matrix & relevant quantities for experiment
    // Calculate nonlinearity parameters
    double g_11 = 4*pi*hbar*hbar/m_1*a_11;
    double g_22 = 4*pi*hbar*hbar/m_2*a_22;
    double g_12 = 2*pi*hbar*hbar*a_12*(m_1+m_2)/(m_1*m_2);
    double g_21 = 2*pi*hbar*hbar*a_21*(m_1+m_2)/(m_1*m_2);

    double Gamma_[2][2];
    Gamma_[0][0] = g_11*T*N_1/(L*L*L*hbar);
    Gamma_[0][1] = g_12*T*N_2/(L*L*L*hbar);
    Gamma_[1][0] = g_21*T*N_1/(L*L*L*hbar);
    Gamma_[1][1] = g_22*T*N_2/(L*L*L*hbar);

    double kappa = L/z_R;
    double delta = 2*L*L/(w_0*w_0);

    //double w = sqrt(w_0*w_0*(1+z*z/(z_R*z_R)));

    double beta_1 = T*L/hbar*m_1*g_F;
    double beta_2 = T*L/hbar*m_2*g_F;

    //Resonance frequencies (_1: Rb, _2: K)
    double omega0[1][2];
    omega0[0][0] = 2*pi*c/r_1;
    omega0[0][1] = 2*pi*c/r_2;

    double Gamma = omega0[0][0]*omega0[0][0]*ee*ee/(6*pi*epsilon0*me*c*c*c);
    double complex alpha_1 = 6*pi*epsilon0*c*c*c*Gamma/ (omega0[0][0]*omega0[0][0]*(omega0[0][0]*omega0[0][0]-omega*omega - I*omega*omega*omega*Gamma/(omega0[0][0]*omega0[0][0])));

    double Gamma_2 = omega0[0][1]*omega0[0][1]*ee*ee/(6*pi*epsilon0*me*c*c*c);
    double complex alpha_2 = 6*pi*epsilon0*c*c*c*Gamma_2/ (omega0[0][1]*omega0[0][1]*(omega0[0][1]*omega0[0][1]-omega*omega - I*omega*omega*omega*Gamma/(omega0[0][1]*omega0[0][1])));

    double eta_1 = P*creal(alpha_1)*T/(pi*hbar*c*epsilon0*w_0*w_0);
    double eta_2 = P*creal(alpha_2)*T/(pi*hbar*c*epsilon0*w_0*w_0);


//---------------------------------------------------------------
//Dipol poteltial: Calculation of offsets and shifts of potential
    printf("\nDipole potential:\n");
    double A = 10.0, B = 2.0;   //Start values
    double A_xml;     //value which will be stored in xml file

    //Potenial with gravity in x direction
    double pot1 = -eta_1*exp(-delta*(x*x+y*y+z*z)/(1+kappa*kappa*z*z))/(1+kappa*kappa*z*z)+beta_1*x; //Fkt.x,y,z,L
    double pot2 = -eta_2*exp(-delta*(x*x+y*y+z*z)/(1+kappa*kappa*z*z))/(1+kappa*kappa*z*z)+beta_2*x; //Fkt.x,y,z,L

    //Determine local minimum
    int n_1 = NewtonMethod_Dipol(A, eta_1, delta, beta_1);
    int n_2 = NewtonMethod_Dipol(B, eta_2, delta, beta_2);
    printf("A\t=\t%E\titeration steps:\t%d\nB\t=\t%E\titeration steps:\t%d\n", A , n_1, B, n_2);
    string output_info_pot;
    double pot_x_plot;  //Plot regime

    //Both species: no local minimum
    if ( ((abs(A)==INFINITY)||(abs(A)>=10E15)||(A==NAN))&&((abs(B)==INFINITY)||(abs(B)>=10E15)||(B==NAN)))
    {
      printf("Attention! no local minimum for BEC_1 and BEC_2: A and B infinite\n");
      output_info_pot ="Attention! No local minimum was found for BEC_1 and BEC_2.";
      pot_x_plot=170;
      A_xml=A;
    }
    //species 2: minimum, species 1: NO local minimum
    else if ((abs(A)==INFINITY)||(abs(A)>=10E15)||(abs(A)==NAN))
    {
      printf("Attention! no local minimum for BEC_1: A infinite\n");
      output_info_pot ="Attention! No local minimum was found for BEC_1.";
      pot_x_plot=(abs(B)+170)*2;
      A_xml=B;
    }
    //species 1: minimum, species 2: NO local minimum
    else if ((abs(B)==INFINITY)||(abs(B)>=10E15)||(abs(B)==NAN))
    {
      printf("Attention! no local minimum for BEC_1: B infinite\n");
      output_info_pot ="Attention! No local minimum was found for BEC_1.";
      pot_x_plot=(abs(A)+170)*2;
      A_xml=A;
    }
    //Both species: local minimum
    else
    {
      output_info_pot ="Local minimum was found.";
      if (abs(A) >= abs(B))
      {
        printf("Local minimum was found!\nShift of potentials: |A| (|A| > |B|)\n");
        A_xml=A;
      }
      else if (abs(B) > abs(A))
      {
        printf("Local minimum was found!\nShift of potentials: |B| (|B| > |A|)\n");
        A_xml=B;
      }
      pot_x_plot=(abs(A)+170)*2;
    }


    //Determine y-Offset:
    double POT_OFFSET = -eta_1*exp(-delta*A*A)+beta_1*A; //Fkt. von L
    double POT_OFFSET2 = -eta_2*exp(-delta*B*B)+beta_2*B; //Fkt. von L
    printf("Pot_Offset:\t%E\t%E\n", POT_OFFSET, POT_OFFSET2);

    //New potentials:
    double pot_1 = -eta_1*exp(-delta*((x+A)*(x+A)+y*y+z*z)/(1+kappa*kappa*z*z))/(1+kappa*kappa*z*z)+beta_1*(x+A)-POT_OFFSET;//Fkt.x,y,z,L
    double pot_2 = -eta_2*exp(-delta*((x+A)*(x+A)+y*y+z*z)/(1+kappa*kappa*z*z))/(1+kappa*kappa*z*z)+beta_2*(x+A)-POT_OFFSET2;//Fkt.x,y,z,L

    //Save new and old Potential in file
    ofstream outf_pot("pot_dipol.txt");
    outf_pot << "#x\tpot_old_1\tpot_new_1\tpot_old_2\tpot_new_2\n";
    int i;
    for (i=-pot_x_plot; i<pot_x_plot; i++)
    {
      pot1 = -eta_1*exp(-delta*i*i)+beta_1*i;
      pot2 = -eta_2*exp(-delta*i*i)+beta_2*i;
      pot_1 = -eta_1*exp(-delta*(i+A)*(i+A))+beta_1*(i+A)-POT_OFFSET;
      pot_2 = -eta_2*exp(-delta*(i+A)*(i+A))+beta_2*(i+A)-POT_OFFSET2;
      outf_pot << i << "\t" << pot1 << "\t" << pot_1 << "\t" << pot2 << "\t" << pot_2 << "\n";
    }


    //---------------------------------------------------
    //Determination of scaling factor of Laplace operator
    double ALPHA_1 = hbar*T/(2*m_1*L*L); //Fkt. von L
    double ALPHA_2 = hbar*T/(2*m_2*L*L); //Fkt. von L



    //---------------------------------------------------------------------------
    //Harmonic trap potental (index ham or ho)
    printf("\nHarmonic trap potental:\n");

    //Frequencie scaling:
    double om_xq_1 = 0.5*T*L*L/hbar*m_1*o_m_x*o_m_x;
    double om_yq_1 = 0.5*T*L*L/hbar*m_1*o_m_y*o_m_y;
    double om_zq_1 = 0.5*T*L*L/hbar*m_1*o_m_z*o_m_z;
    double om_xq_2 = 0.5*T*L*L/hbar*m_2*o_m_x*o_m_x;
    double om_yq_2 = 0.5*T*L*L/hbar*m_2*o_m_y*o_m_y;
    double om_zq_2 = 0.5*T*L*L/hbar*m_2*o_m_z*o_m_z;


    double A_ham=-10.0, B_ham=-10.0;  //start value
    double A_ham_xml;     //value which will be stored in xml file

    //Potenial with gravity in x direction
    double pot_ham1 = om_xq_1*x*x+om_yq_1*y*y+om_zq_1*z*z+beta_1*x; //Fkt. von x,y,z
    double pot_ham2 = om_xq_2*x*x+om_yq_2*y*y+om_zq_2*z*z+beta_2*x; //Fkt. von x,y,z

    int n_ham_1=NewtonMethod_ham_osci(A_ham, om_xq_1, beta_1);
    int n_ham_2=NewtonMethod_ham_osci(B_ham, om_xq_2, beta_2);
    printf("A_ho\t=\t%E\titeration steps:\t%d\nB_ho\t=\t%E\titeration steps:\t%d\n", A_ham , n_ham_1, B_ham, n_ham_2);
    string output_info_pot_ham;
    double pot_x_plot_ham;  //Plot regime

    if ( ((abs(A_ham)==INFINITY)||(abs(A_ham)>=10E15)||(A_ham==NAN))&&((abs(B_ham)==INFINITY)||(abs(B_ham)>=10E15)||(B_ham==NAN)))
    {
      printf("Attention! no local minimum for BEC_1 and BEC_2: A_ho and B_ho infinite\n");
      output_info_pot ="Attention! No local minimum was found for BEC_1 and BEC_2.";
      pot_x_plot_ham= 20;
      A_ham_xml=A_ham;
    }
    else if ((abs(A_ham)==INFINITY)||(abs(A_ham)>=10E15)||(abs(A_ham)==NAN))
    {
      printf("Attention! no local minimum for BEC_1: A_ho infinite\n");
      output_info_pot ="Attention! No local minimum was found for BEC_1.";
      pot_x_plot_ham= (abs(B_ham)+20)*2;
      A_ham_xml=B_ham;
    }
    else if ((abs(B_ham)==INFINITY)||(abs(B_ham)>=10E15)||(abs(B_ham)==NAN))
    {
      printf("Attention! no local minimum for BEC_1: B_ho infinite\n");
      output_info_pot ="Attention! No local minimum was found for BEC_1.";
      pot_x_plot_ham= (abs(A_ham)+20)*2;
      A_ham_xml=A_ham;
    }
    else
    {
      output_info_pot_ham ="Local minimum was found";
      if (abs(A_ham) >= abs(B_ham))
      {
        printf("Local minimum was found!\nShift of potentials: |A_ho| (|A_ho| > |B_ho|)\n");
        A_ham_xml=A_ham;
      }
      else if (abs(B_ham) > abs(A_ham))
      {
        printf("Local minimum was found!\nShift of potentials: |B_ho| (|B_ho| > |A_ho|)\n");
        A_ham_xml=B_ham;
      }
      pot_x_plot_ham=(abs(A_ham)+20)*2;
    }

    //Determine y-Offset
    double POT_OFFSET_ham = om_xq_1*A_ham*A_ham+beta_1*A_ham;
    double POT_OFFSET_ham2 = om_xq_2*B_ham*B_ham+beta_2*B_ham;
    printf("Pot_Offset:\t%E\t%E\n", POT_OFFSET_ham, POT_OFFSET_ham2);

    //Save new and old Potential in file
    ofstream outf_pot_ham("pot_ham.txt");
    outf_pot_ham << "#x\tpot_old_1\tpot_new_1\tpot_old_2\tpot_new_2\n";

    double pot_ham_1;
    double pot_ham_2;
    //for(i=-abs(A_ham)*2; i<abs(A_ham)*2; i++)

    for (i=-pot_x_plot_ham; i<pot_x_plot_ham; i++)
    {
      pot_ham_1 = om_xq_1*(i+A_ham)*(i+A_ham)+beta_1*(i+A_ham)-POT_OFFSET_ham; //Fkt. von x
      pot_ham_2 = om_xq_2*(i+A_ham)*(i+A_ham)+beta_2*(i+A_ham)-POT_OFFSET_ham2; //Fkt. von x
      pot_ham1 = om_xq_1*i*i+beta_1*i;
      pot_ham2 = om_xq_2*i*i+beta_2*i;
      outf_pot_ham << i << "\t" << pot_ham1 << "\t" << pot_ham_1 << "\t" << pot_ham2 << "\t" << pot_ham_2 << "\n";
    }

    //----------------
    //Save in xml-file
    if (species==12)
    {
      printf("\nTwo Species\n");
      ofstream outf( "BEC_12.xml" );
      outf << setprecision(10);
      outf << "<SET_PARAMS>" << "\n";
      outf << "<!-- T="<< T <<", L="<< L <<", N1 = "<< N_1 <<", N2= "<< N_2<<", P="<< P <<", omega_x="<< o_m_x <<", omega_y="<< o_m_y <<", omega_z="<< o_m_z << ", w0="<< w_0 <<", lambda="<< lambda <<"-->\n";
      outf << "<DATEI_NAME>" << filename << "</DATEI_NAME>"<< "\n";
      outf << "<DATEI_NAME_2>" << filename2 << "</DATEI_NAME_2>"<< "\n";
      outf << "<FIN_N val=\"" << N_1/N_1 << "\" val_2=\"" << N_2/N_2 << "\"></FIN_N>"<< "\n";
      outf << "<MXITER val=\"" << MaxIter << "\"></MXITER>"<< "\n";
      outf << "<STEP_SIZE val=\"" << stepsize << "\"></STEP_SIZE>"<< "\n";
      outf << "<EPS val=\"" << eps << "\"></EPS>"<< "\n";
      outf << "<NA val=\"" << Na << "\"></NA>"<< "\n";
      outf << "<NK val=\"" << Nk << "\"></NK>"<< "\n";
      outf << "<GS val=\"" << Gamma_[0][0] << "\" val_12=\"" << Gamma_[0][1] << "\" val_21=\"" << Gamma_[1][0] << "\" val_2=\"" << Gamma_[1][1] << "\"></GS>"<< "\n";
      outf << "<ALPHA val=\"" << ALPHA_1 << "\" val_2=\"" << ALPHA_2 << "\"></ALPHA>"<< "\n";
      outf << "<BETA val=\"" << beta_1 << "\" val_2=\"" << beta_2 << "\"></BETA>"<< "\n";
      outf << "<ETA val=\"" << eta_1 << "\" val_2=\"" << eta_2 << "\"></ETA>"<< "\n";
      outf << "<KAPPA val=\"" << kappa << "\"></KAPPA>"<< "\n";
      outf << "<DELTA val=\"" << delta << "\"></DELTA>"<< "\n";
      outf << "<T_OFF val=\"" << t_off << "\"></T_OFF>"<< "\n";
      outf << "<T_ON val=\"" << t_on << "\"></T_ON>"<< "\n";
      outf << "<!-- Dipol potential: "<< output_info_pot << "-->" << "\n";
      outf << "<A val=\"" << A_xml << "\"></A>"<< "<!-- x shift -->\n";
      outf << "<POT_OFFSET val=\"" << POT_OFFSET << "\" val_2=\"" << POT_OFFSET2 << "\"></POT_OFFSET>"<< "\n";
      outf << "<!-- Harmonic oscillator potential: "<< output_info_pot_ham << "-->" << "\n";
      outf << "<FF om_xq=\"" << om_xq_1 << "\" om_yq=\"" << om_yq_1 << "\" om_zq=\"" << om_zq_1 << "\"></FF>"<< "\n";
      outf << "<FF2 om_xq_2=\"" << om_xq_2 << "\" om_yq_2=\"" << om_yq_2 << "\" om_zq_2=\"" << om_zq_2 << "\"></FF2>"<< "\n";
      outf << "<B val=\"" << A_ham_xml << "\"></B>"<< "<!-- x shift -->\n";
      outf << "<POT_OFFSET_HO val=\"" << POT_OFFSET_ham << "\" val_2=\"" << POT_OFFSET_ham2 << "\"></POT_OFFSET_HO>"<< "\n";
      outf << "</SET_PARAMS>" << "\n";
    }
    else if (species==1)
    {
      printf("\nFirst Species\n");
      ofstream outf( "BEC_1.xml" );
      outf << setprecision(10);
      outf << "<SET_PARAMS>" << "\n";
      outf << "<!-- T="<< T <<", L="<< L <<", N1 = "<< N_1 <<", N2= "<< N_2<<", P="<< P <<", omega_r_1="<< o_m_x <<", omega_r_2="<< o_m_y <<", omega_z_1="<< o_m_z <<", omega_z_2="<< o_m_z <<", w0="<< w_0 <<", lambda="<< lambda <<"-->\n";
      outf << "<DATEI_NAME>" << filename << "</DATEI_NAME>"<< "\n";
      outf << "<FIN_N val=\"" << N_1/N_1 << "\"></FIN_N>"<< "\n";
      outf << "<MXITER val=\"" << MaxIter << "\"></MXITER>"<< "\n";
      outf << "<STEP_SIZE val=\"" << stepsize << "\"></STEP_SIZE>"<< "\n";
      outf << "<EPS val=\"" << eps << "\"></EPS>"<< "\n";
      outf << "<NA val=\"" << Na << "\"></NA>"<< "\n";
      outf << "<NK val=\"" << Nk << "\"></NK>"<< "\n";
      outf << "<GS val=\"" << Gamma_[0][0] << "\"></GS>"<< "\n";
      outf << "<ALPHA val=\"" << ALPHA_1 << "\"></ALPHA>"<< "\n";
      outf << "<BETA val=\"" << beta_1 << "\"></BETA>"<< "\n";
      outf << "<ETA val=\"" << eta_1 << "\"></ETA>"<< "\n";
      outf << "<KAPPA val=\"" << kappa << "\"></KAPPA>"<< "\n";
      outf << "<DELTA val=\"" << delta << "\"></DELTA>"<< "\n";
      outf << "<T_OFF val=\"" << t_off << "\"></T_OFF>"<< "\n";
      outf << "<T_ON val=\"" << t_on << "\"></T_ON>"<< "\n";
      outf << "<!-- Dipol potential: "<< output_info_pot << "-->" << "\n";
      outf << "<A val=\"" << A << "\"></A>"<< "<!-- x shift -->\n";
      outf << "<POT_OFFSET val=\"" << POT_OFFSET << "\"></POT_OFFSET>"<< "\n";
      outf << "<!-- Harmonic oscillator potental: "<< output_info_pot_ham << "-->" << "\n";
      outf << "<FF om_xq=\"" << om_xq_1 << "\" om_yq=\"" << om_yq_1 << "\" om_zq=\"" << om_zq_1 << "\"></FF>"<< "\n";
      outf << "<B val=\"" << A_ham << "\"></B>"<< "<!-- x shift -->\n";
      outf << "<POT_OFFSET_HO val=\"" << POT_OFFSET_ham << "\"></POT_OFFSET_HO>"<< "\n";
      outf << "</SET_PARAMS>" << "\n";
    }
    else if (species==2)
    {
      printf("\nSecond Species\n");
      ofstream outf( "BEC_2.xml" );
      outf << setprecision(10);
      outf << "<SET_PARAMS>" << "\n";
      outf << "<!-- T="<< T <<", L="<< L <<", N1 = "<< N_1 <<", N2= "<< N_2<<", P="<< P <<", omega_r_1="<< o_m_x <<", omega_r_2="<< o_m_y <<", omega_z_1="<< o_m_z <<", omega_z_2="<< o_m_z <<", w0="<< w_0 <<", lambda="<< lambda <<"-->\n";
      outf << "<DATEI_NAME>" << filename2 << "</DATEI_NAME>"<< "\n";
      outf << "<FIN_N val=\"" << N_2/N_2 << "\"></FIN_N>"<< "\n";
      outf << "<MXITER val=\"" << MaxIter << "\"></MXITER>"<< "\n";
      outf << "<STEP_SIZE val=\"" << stepsize << "\"></STEP_SIZE>"<< "\n";
      outf << "<EPS val=\"" << eps << "\"></EPS>"<< "\n";
      outf << "<NA val=\"" << Na << "\"></NA>"<< "\n";
      outf << "<NK val=\"" << Nk << "\"></NK>"<< "\n";
      outf << "<GS val=\"" << Gamma_[1][1] << "\"></GS>"<< "\n";
      outf << "<ALPHA val=\"" << ALPHA_2 << "\"></ALPHA>"<< "\n";
      outf << "<BETA val=\"" << beta_2 << "\"></BETA>"<< "\n";
      outf << "<ETA val=\"" << eta_2 << "\"></ETA>"<< "\n";
      outf << "<KAPPA val=\"" << kappa << "\"></KAPPA>"<< "\n";
      outf << "<DELTA val=\"" << delta << "\"></DELTA>"<< "\n";
      outf << "<T_OFF val=\"" << t_off << "\"></T_OFF>"<< "\n";
      outf << "<T_ON val=\"" << t_on << "\"></T_ON>"<< "\n";
      outf << "<!-- Dipol potential: "<< output_info_pot << "-->" << "\n";
      outf << "<A val=\"" << B << "\"></A>"<< "<!-- x shift -->\n";
      outf << "<POT_OFFSET val=\"" << POT_OFFSET2 << "\"></POT_OFFSET>"<< "\n";
      outf << "<!-- Harmonic oscillator potential: "<< output_info_pot_ham << "-->" << "\n";
      outf << "<FF om_xq_2=\"" << om_xq_2 << "\" om_yq_2=\"" << om_yq_2 << "\" om_zq_2=\"" << om_zq_2 << "\"></FF>"<< "\n";
      outf << "<B val=\"" << B_ham << "\"></B>"<< "<!-- x shift -->\n";
      outf << "<POT_OFFSET_HO val=\"" << POT_OFFSET_ham2 << "\"></POT_OFFSET_HO>"<< "\n";
      outf << "</SET_PARAMS>" << "\n";
    }
  }
}
