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
#include <cmath>
#include <omp.h>
#include "cft_2d_MPI.h"
#include "cft_3d_MPI.h"
#include "muParser.h"
#include "ParameterHandler.h"
#include "CRT_Base_IF_2_mpi.h"

using namespace std;

namespace MPI
{
  namespace RT_Solver_2
  {
    template<class T,int dim>
    class Bragg_double : public CRT_Base_IF_2_mpi<T,dim,6>
    {
    public:
      Bragg_double( ParameterHandler * );
      virtual ~Bragg_double();

    protected:
      void Do_Bragg_ad();
      void Save_Psi_xy();

      static void Do_Bragg_ad_Wrapper(void *,sequence_item &);
      static void Save_Psi_xy_Wrapper(void *,sequence_item &);

      bool run_custom_sequence( const sequence_item & );
      bool check_consistency( const generic_header &head );

      double *m_Mirror;
      int m_no_of_Mirror_pts;

      using CRT_Base_2_mpi<T,dim,6>::m_map_stepfcts;
      using CRT_Base_2_mpi<T,dim,6>::m_fields;
      using CRT_Base_2_mpi<T,dim,6>::m_header;
      using CRT_Base_2_mpi<T,dim,6>::m_gs;
      using CRT_Base_2_mpi<T,dim,6>::MTime;
      using CRT_Base_IF_2_mpi<T,dim,6>::Amp;
      using CRT_Base_IF_2_mpi<T,dim,6>::Amp2;
      using CRT_Base_IF_2_mpi<T,dim,6>::laser_k;
      using CRT_Base_IF_2_mpi<T,dim,6>::laser_dk;
      using CRT_Base_IF_2_mpi<T,dim,6>::laser_k2;
      using CRT_Base_IF_2_mpi<T,dim,6>::laser_dk2;
      using CRT_Base_IF_2_mpi<T,dim,6>::chirp_rate;
      using CRT_Base_IF_2_mpi<T,dim,6>::laser_domh;
      using CRT_Base_IF_2_mpi<T,dim,6>::laser_domh2;
      using CRT_Base_IF_2_mpi<T,dim,6>::beta;
      using CRT_Base_IF_2_mpi<T,dim,6>::phase;
      using CRT_Base_IF_2_mpi<T,dim,6>::DeltaL;
    };

    template<class T,int dim>
    Bragg_double<T,dim>::Bragg_double( ParameterHandler *p ) : CRT_Base_IF_2_mpi<T,dim,6>( p )
    {
      m_map_stepfcts["bragg_ad"] = &Do_Bragg_ad_Wrapper;

      CPoint<dim> pt1;
      CPoint<dim> pt2;
      pt2[0] = 2*this->laser_k[0];
      CPoint<dim> pt3;
      pt3[0] = -2*this->laser_k[0];
      this->m_rabi_momentum_list.push_back(pt1);
      this->m_rabi_momentum_list.push_back(pt2);
      this->m_rabi_momentum_list.push_back(pt3);

      pt1[0] = 0;
      pt2[0] = 2*this->laser_k2[0];
      pt3[0] = -2*this->laser_k2[0];

      this->m_rabi_momentum_list2.push_back(pt1);
      this->m_rabi_momentum_list2.push_back(pt2);
      this->m_rabi_momentum_list2.push_back(pt3);

      switch (dim)
      {
      case 1:
        m_no_of_Mirror_pts = this->Get_dimY();
        break;
      case 2:
        m_no_of_Mirror_pts = this->Get_dimY()*this->Get_dimZ();
        this->m_custom_fct=&Save_Psi_xy_Wrapper;
        break;
      }

      m_Mirror = new double[m_no_of_Mirror_pts];
      memset( m_Mirror, 0, sizeof(double)*m_no_of_Mirror_pts);
    }

    template<class T,int dim>
    Bragg_double<T,dim>::~Bragg_double()
    {
      delete [] m_Mirror;
    }

    template<class T,int dim>
    bool Bragg_double<T,dim>::check_consistency( const generic_header &head )
    {
      bool retval=false;
      switch ( dim )
      {
      case 2:
        retval = ( head.nDimX == this->Get_dimY() );
        break;
      case 3:
        retval = ( head.nDimX == this->Get_dimY() ) && ( head.nDimY == this->Get_dimZ() );
        break;
      }
      return retval;
    }

    template<class T,int dim>
    bool Bragg_double<T,dim>::run_custom_sequence( const sequence_item &item )
    {
      if ( item.name == "set_mirror" )
      {
        if ( this->m_myrank == 0 )
        {
          std::cout << "Mirror" <<std::endl;
          cout << item.content << endl;
          generic_header header_tmp;
          ifstream file1(item.content, ifstream::binary );
          file1.read( (char *)&header_tmp, sizeof(generic_header));

          if ( check_consistency(header_tmp) )
          {
            file1.seekg( sizeof(generic_header), ifstream::beg );
            file1.read( (char *)m_Mirror, sizeof(double)*m_no_of_Mirror_pts );
            file1.close();
          }
          else
          {
            cout << "Change number of points for mirror." << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
          }
        }
        MPI_Bcast( m_Mirror, m_no_of_Mirror_pts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        return true;
      }
      return false;
    }

    template<class T,int dim>
    void Bragg_double<T,dim>::Save_Psi_xy_Wrapper ( void *ptr, sequence_item & /*item*/ )
    {
      Bragg_double *self = static_cast<Bragg_double *>(ptr);
      self->Save_Psi_xy();
    }

    template<class T,int dim>
    void Bragg_double<T,dim>::Save_Psi_xy()
    {
      const int comp=0;
      char filename[1024];

      ptrdiff_t dimY = this->Get_dimY();
      ptrdiff_t dimZ = this->Get_dimZ();
      ptrdiff_t k= dimZ/2;

      fftw_complex *Psi_tmp = fftw_alloc_complex( this->m_loc_dimX*dimY );
      fftw_complex *Psi = m_fields[comp]->Get_p2_Data();

      for ( ptrdiff_t i=0; i<this->m_loc_dimX; i++ )
      {
        for ( ptrdiff_t j=0; j<dimY; j++ )
        {
          Psi_tmp[j+dimY*i][0] = Psi[k+dimZ*(j+dimY*i)][0];
          Psi_tmp[j+dimY*i][1] = Psi[k+dimZ*(j+dimY*i)][1];
        }
      }

      double *data = reinterpret_cast<double *>(Psi_tmp);

      MPI_Status status;
      MPI_File   fh;

      m_header.nDims = 2;
      m_header.nDimZ = 1;

      MPI_Offset offset = sizeof(generic_header) + sizeof(fftw_complex)*this->m_loc_start_dimX*dimY;

      sprintf( filename, "S3d_%.3f_%d.bin", this->Get_t(), comp+1 );
      MPI_File_open( MPI_COMM_WORLD, const_cast<char *>(filename), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh );

      if ( this->m_myrank == 0 )
        MPI_File_write( fh, &m_header, sizeof(generic_header), MPI_BYTE, &status );

      MPI_File_write_at( fh, offset, data, 2*this->m_loc_dimX*dimY, MPI_DOUBLE, MPI_STATUS_IGNORE );
      MPI_File_close( &fh );

      m_header.nDims = 3;
      m_header.nDimZ = dimZ;
    }

    template<class T,int dim>
    void Bragg_double<T,dim>::Do_Bragg_ad_Wrapper ( void *ptr, sequence_item & /*item*/ )
    {
      Bragg_double *self = static_cast<Bragg_double *>(ptr);
      self->Do_Bragg_ad();
    }

    template<class T,int dim>
    void Bragg_double<T,dim>::Do_Bragg_ad()
    {
      //MTime.enter_section("Do_NL_Step_light_ana");

      // first species

      fftw_complex *Psi_1 = this->m_fields[0]->Get_p2_Data();
      fftw_complex *Psi_2 = this->m_fields[1]->Get_p2_Data();
      fftw_complex *Psi_3 = this->m_fields[2]->Get_p2_Data();
      fftw_complex *Psi_4 = this->m_fields[3]->Get_p2_Data();
      fftw_complex *Psi_5 = this->m_fields[4]->Get_p2_Data();
      fftw_complex *Psi_6 = this->m_fields[5]->Get_p2_Data();
      T *ft = reinterpret_cast<T *>(this->m_fields[0]);

      const double dt = this->Get_dt();
      const double t1 = this->Get_t()+0.5*dt;
      CPoint<dim> x;
      vector<ptrdiff_t> global_idx(3);

      fftw_complex O11, O12, O21, O22, O13, O31, O33, O32, O23, gamma_1, gamma_2, gamma_3, eta;
      double re1, im1, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, V11, V22, E1, E2, Omega_p, Omega_m;
      ptrdiff_t jk;

      for ( ptrdiff_t l=0; l<this->m_no_of_pts; l++ )
      {
        x = ft->Get_x(l);

        switch ( dim )
        {
        case 2:
          ft->local_to_global_rs( l, global_idx );
          jk = global_idx[1];
          break;
        case 3:
          ft->local_to_global_rs( l, global_idx );
          jk = global_idx[2]+this->Get_dimY()*global_idx[1];
          break;
        }

        tmp1 = Psi_1[l][0]*Psi_1[l][0]+Psi_1[l][1]*Psi_1[l][1];
        tmp2 = Psi_2[l][0]*Psi_2[l][0]+Psi_2[l][1]*Psi_2[l][1];
        tmp3 = Psi_3[l][0]*Psi_3[l][0]+Psi_3[l][1]*Psi_3[l][1];
        tmp4 = Psi_4[l][0]*Psi_4[l][0]+Psi_4[l][1]*Psi_4[l][1];
        tmp5 = Psi_5[l][0]*Psi_5[l][0]+Psi_5[l][1]*Psi_5[l][1];
        tmp6 = Psi_6[l][0]*Psi_6[l][0]+Psi_6[l][1]*Psi_6[l][1];

        V11 = m_gs[0]*tmp1 + m_gs[1]*tmp2 + m_gs[2]*tmp3 + m_gs[3]*tmp4 + m_gs[4]*tmp5 + m_gs[5]*tmp6 + beta*x;
        V22 = m_gs[6]*tmp1 + m_gs[7]*tmp2 + m_gs[8]*tmp3 + m_gs[9]*tmp4 + m_gs[10]*tmp5 + m_gs[11]*tmp6 - DeltaL[1]+beta*x;

        Omega_p = Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0]+chirp_rate[0]*t1)*t1);
        Omega_m = Amp[0]*cos(-laser_k[0]*x[0]-(laser_domh[0]-chirp_rate[0]*t1)*t1);
        //Omega_p = Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0]+chirp_rate[0]*t1)*t1-0.5*(phase_L[0]+m_Mirror[jk]));
        //Omega_m = Amp[0]*cos(-laser_k[0]*x[0]-(laser_domh[0]-chirp_rate[0]*t1)*t1-0.5*(phase_L[0]-m_Mirror[jk]));
        //Omega_p = Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0])*t1);
        //Omega_p = Amp[0]*cos(-laser_k[0]*x[0]-(laser_domh[0])*t1);

        if ((Omega_m == 0.0 ) && ( Omega_p == 0.0 ))
        {
          sincos( -dt*V11, &im1, &re1 );
          tmp1 = Psi_1[l][0];
          Psi_1[l][0] = tmp1*re1-Psi_1[l][1]*im1;
          Psi_1[l][1] = tmp1*im1+Psi_1[l][1]*re1;

          sincos( -dt*V22, &im1, &re1 );
          tmp1 = Psi_2[l][0];
          Psi_2[l][0] = tmp1*re1-Psi_2[l][1]*im1;
          Psi_2[l][1] = tmp1*im1+Psi_2[l][1]*re1;
          tmp1 = Psi_3[l][0];
          Psi_3[l][0] = tmp1*re1-Psi_3[l][1]*im1;
          Psi_3[l][1] = tmp1*im1+Psi_3[l][1]*re1;
          continue;
        }

        sincos( -0.5*laser_dk[0]*x[0], &im1, &re1 );
        eta[0] = re1;
        eta[1] = im1;

        tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*(Omega_p*Omega_p+Omega_m*Omega_m));
        E1 = 0.5*(V11+V22+tmp1);
        E2 = 0.5*(V11+V22-tmp1);

        sincos( -dt*E1, &im1, &re1 );
        tmp1 = 1.0/fabs((V22-E1)*(V22-E1)+Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_1[0] = tmp1*re1;
        gamma_1[1] = tmp1*im1;

        sincos( -dt*E2, &im1, &re1 );
        tmp1 = 1.0/fabs((V22-E2)*(V22-E2)+Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_2[0] = tmp1*re1;
        gamma_2[1] = tmp1*im1;

        sincos( -dt*V22, &im1, &re1);
        tmp1 = 1.0/fabs(Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_3[0] = tmp1*re1;
        gamma_3[1] = tmp1*im1;

        tmp1 = (V22-E2)*(V22-E2);
        tmp2 = (V22-E1)*(V22-E1);
        O11[0] = tmp1*gamma_2[0] + tmp2*gamma_1[0];
        O11[1] = tmp1*gamma_2[1] + tmp2*gamma_1[1];

        tmp1 = (V22-E2)*gamma_2[0] + (V22-E1)*gamma_1[0];
        tmp2 = (V22-E2)*gamma_2[1] + (V22-E1)*gamma_1[1];

        O12[0] = -Omega_p*tmp1;
        O12[1] = -Omega_p*tmp2;

        O13[0] = -Omega_m*tmp1;
        O13[1] = -Omega_m*tmp2;

        O21[0] = O12[0];
        O21[1] = O12[1];

        O31[0] = O13[0];
        O31[1] = O13[1];

        tmp1 = O12[0];
        O12[0] = eta[0]*tmp1-eta[1]*O12[1];
        O12[1] = eta[1]*tmp1+eta[0]*O12[1];

        tmp1 = O21[0];
        O21[0] = eta[0]*tmp1+eta[1]*O21[1];
        O21[1] = eta[0]*O21[1]-eta[1]*tmp1;

        tmp1 = O31[0];
        O31[0] = eta[0]*tmp1-eta[1]*O31[1];
        O31[1] = eta[1]*tmp1+eta[0]*O31[1];

        tmp1 = O13[0];
        O13[0] = eta[0]*tmp1+eta[1]*O13[1];
        O13[1] = eta[0]*O13[1]-eta[1]*tmp1;

        tmp1 = Omega_p*Omega_p;
        tmp2 = Omega_m*Omega_m;
        O22[0] = tmp1*(gamma_2[0]+gamma_1[0]) + tmp2*gamma_3[0];
        O22[1] = tmp1*(gamma_2[1]+gamma_1[1]) + tmp2*gamma_3[1];

        O33[0] = tmp2*(gamma_2[0]+gamma_1[0]) + tmp1*gamma_3[0];
        O33[1] = tmp2*(gamma_2[1]+gamma_1[1]) + tmp1*gamma_3[1];

        tmp1 = Omega_m*Omega_p;
        O23[0] = tmp1*(gamma_1[0]+gamma_2[0]-gamma_3[0]);
        O23[1] = tmp1*(gamma_1[1]+gamma_2[1]-gamma_3[1]);

        O32[0] = O23[0];
        O32[1] = O23[1];

        tmp1 = O23[0];
        O23[0] = (eta[0]*eta[0]+eta[1]*eta[1])*tmp1+2*eta[0]*eta[1]*O23[1];
        O23[1] = (eta[0]*eta[0]+eta[1]*eta[1])*O23[1]-2*eta[0]*eta[1]*tmp1;

        tmp1 = O32[0];
        O32[0] = (eta[0]*eta[0]-eta[1]*eta[1])*tmp1-2*eta[0]*eta[1]*O32[1];
        O32[1] = (eta[0]*eta[0]-eta[1]*eta[1])*O23[1]+2*eta[0]*eta[1]*tmp1;

        gamma_1[0] = Psi_1[l][0];
        gamma_1[1] = Psi_1[l][1];
        gamma_2[0] = Psi_2[l][0];
        gamma_2[1] = Psi_2[l][1];
        gamma_3[0] = Psi_3[l][0];
        gamma_3[1] = Psi_3[l][1];

        Psi_1[l][0] = (O11[0]*gamma_1[0]-O11[1]*gamma_1[1]) + (O12[0]*gamma_2[0]-O12[1]*gamma_2[1]) + (O13[0]*gamma_3[0]-O13[1]*gamma_3[1]);
        Psi_1[l][1] = (O11[0]*gamma_1[1]+O11[1]*gamma_1[0]) + (O12[0]*gamma_2[1]+O12[1]*gamma_2[0]) + (O13[0]*gamma_3[1]+O13[1]*gamma_3[0]);

        Psi_2[l][0] = (O21[0]*gamma_1[0]-O21[1]*gamma_1[1]) + (O22[0]*gamma_2[0]-O22[1]*gamma_2[1]) + (O23[0]*gamma_3[0]-O23[1]*gamma_3[1]);
        Psi_2[l][1] = (O21[0]*gamma_1[1]+O21[1]*gamma_1[0]) + (O22[0]*gamma_2[1]+O22[1]*gamma_2[0]) + (O23[0]*gamma_3[1]+O23[1]*gamma_3[0]);

        Psi_3[l][0] = (O31[0]*gamma_1[0]-O31[1]*gamma_1[1]) + (O32[0]*gamma_2[0]-O32[1]*gamma_2[1]) + (O33[0]*gamma_3[0]-O33[1]*gamma_3[1]);
        Psi_3[l][1] = (O31[0]*gamma_1[1]+O31[1]*gamma_1[0]) + (O32[0]*gamma_2[1]+O32[1]*gamma_2[0]) + (O33[0]*gamma_3[1]+O33[1]*gamma_3[0]);
      }

      // second species
      Psi_1 = this->m_fields[3]->Get_p2_Data();
      Psi_2 = this->m_fields[4]->Get_p2_Data();
      Psi_3 = this->m_fields[5]->Get_p2_Data();
      Psi_4 = this->m_fields[0]->Get_p2_Data();
      Psi_5 = this->m_fields[1]->Get_p2_Data();
      Psi_6 = this->m_fields[2]->Get_p2_Data();

      for ( ptrdiff_t l=0; l<this->m_no_of_pts; l++ )
      {
        x = ft->Get_x(l);

        switch ( dim )
        {
        case 2:
          ft->local_to_global_rs( l, global_idx );
          jk = global_idx[1];
          break;
        case 3:
          ft->local_to_global_rs( l, global_idx );
          jk = global_idx[2]+this->Get_dimY()*global_idx[1];
          break;
        }

        tmp1 = Psi_1[l][0]*Psi_1[l][0]+Psi_1[l][1]*Psi_1[l][1];
        tmp2 = Psi_2[l][0]*Psi_2[l][0]+Psi_2[l][1]*Psi_2[l][1];
        tmp3 = Psi_3[l][0]*Psi_3[l][0]+Psi_3[l][1]*Psi_3[l][1];
        tmp4 = Psi_4[l][0]*Psi_4[l][0]+Psi_4[l][1]*Psi_4[l][1];
        tmp5 = Psi_5[l][0]*Psi_5[l][0]+Psi_5[l][1]*Psi_5[l][1];
        tmp6 = Psi_6[l][0]*Psi_6[l][0]+Psi_6[l][1]*Psi_6[l][1];

        V11 = m_gs[43]*tmp1 + m_gs[44]*tmp2 + m_gs[45]*tmp3 + m_gs[3]*tmp4 + m_gs[4]*tmp5 + m_gs[5]*tmp6 + beta*x;
        V22 = m_gs[53]*tmp1 + m_gs[54]*tmp2 + m_gs[55]*tmp3 + m_gs[9]*tmp4 + m_gs[10]*tmp5 + m_gs[11]*tmp6 - DeltaL[1]+beta*x;

        Omega_p = Amp2[0]*cos(laser_k2[0]*x[0]-(laser_domh2[0]+chirp_rate[0]*t1)*t1);
        Omega_m = Amp2[0]*cos(-laser_k2[0]*x[0]-(laser_domh2[0]-chirp_rate[0]*t1)*t1);
        //Omega_p = Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0]+chirp_rate[0]*t1)*t1-0.5*(phase_L[0]+m_Mirror[jk]));
        //Omega_m = Amp[0]*cos(-laser_k[0]*x[0]-(laser_domh[0]-chirp_rate[0]*t1)*t1-0.5*(phase_L[0]-m_Mirror[jk]));
        //Omega_p = Amp[0]*cos(laser_k[0]*x[0]-(laser_domh[0])*t1);
        //Omega_p = Amp[0]*cos(-laser_k[0]*x[0]-(laser_domh[0])*t1);

        if ((Omega_m == 0.0 ) && ( Omega_p == 0.0 ))
        {
          sincos( -dt*V11, &im1, &re1 );
          tmp1 = Psi_1[l][0];
          Psi_1[l][0] = tmp1*re1-Psi_1[l][1]*im1;
          Psi_1[l][1] = tmp1*im1+Psi_1[l][1]*re1;

          sincos( -dt*V22, &im1, &re1 );
          tmp1 = Psi_2[l][0];
          Psi_2[l][0] = tmp1*re1-Psi_2[l][1]*im1;
          Psi_2[l][1] = tmp1*im1+Psi_2[l][1]*re1;
          tmp1 = Psi_3[l][0];
          Psi_3[l][0] = tmp1*re1-Psi_3[l][1]*im1;
          Psi_3[l][1] = tmp1*im1+Psi_3[l][1]*re1;
          continue;
        }

        sincos( -0.5*laser_dk[0]*x[0], &im1, &re1 );
        eta[0] = re1;
        eta[1] = im1;

        tmp1 = sqrt((V11-V22)*(V11-V22)+4.0*(Omega_p*Omega_p+Omega_m*Omega_m));
        E1 = 0.5*(V11+V22+tmp1);
        E2 = 0.5*(V11+V22-tmp1);

        sincos( -dt*E1, &im1, &re1 );
        tmp1 = 1.0/fabs((V22-E1)*(V22-E1)+Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_1[0] = tmp1*re1;
        gamma_1[1] = tmp1*im1;

        sincos( -dt*E2, &im1, &re1 );
        tmp1 = 1.0/fabs((V22-E2)*(V22-E2)+Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_2[0] = tmp1*re1;
        gamma_2[1] = tmp1*im1;

        sincos( -dt*V22, &im1, &re1);
        tmp1 = 1.0/fabs(Omega_p*Omega_p+Omega_m*Omega_m);
        gamma_3[0] = tmp1*re1;
        gamma_3[1] = tmp1*im1;

        tmp1 = (V22-E2)*(V22-E2);
        tmp2 = (V22-E1)*(V22-E1);
        O11[0] = tmp1*gamma_2[0] + tmp2*gamma_1[0];
        O11[1] = tmp1*gamma_2[1] + tmp2*gamma_1[1];

        tmp1 = (V22-E2)*gamma_2[0] + (V22-E1)*gamma_1[0];
        tmp2 = (V22-E2)*gamma_2[1] + (V22-E1)*gamma_1[1];

        O12[0] = -Omega_p*tmp1;
        O12[1] = -Omega_p*tmp2;

        O13[0] = -Omega_m*tmp1;
        O13[1] = -Omega_m*tmp2;

        O21[0] = O12[0];
        O21[1] = O12[1];

        O31[0] = O13[0];
        O31[1] = O13[1];

        tmp1 = O12[0];
        O12[0] = eta[0]*tmp1-eta[1]*O12[1];
        O12[1] = eta[1]*tmp1+eta[0]*O12[1];

        tmp1 = O21[0];
        O21[0] = eta[0]*tmp1+eta[1]*O21[1];
        O21[1] = eta[0]*O21[1]-eta[1]*tmp1;

        tmp1 = O31[0];
        O31[0] = eta[0]*tmp1-eta[1]*O31[1];
        O31[1] = eta[1]*tmp1+eta[0]*O31[1];

        tmp1 = O13[0];
        O13[0] = eta[0]*tmp1+eta[1]*O13[1];
        O13[1] = eta[0]*O13[1]-eta[1]*tmp1;

        tmp1 = Omega_p*Omega_p;
        tmp2 = Omega_m*Omega_m;
        O22[0] = tmp1*(gamma_2[0]+gamma_1[0]) + tmp2*gamma_3[0];
        O22[1] = tmp1*(gamma_2[1]+gamma_1[1]) + tmp2*gamma_3[1];

        O33[0] = tmp2*(gamma_2[0]+gamma_1[0]) + tmp1*gamma_3[0];
        O33[1] = tmp2*(gamma_2[1]+gamma_1[1]) + tmp1*gamma_3[1];

        tmp1 = Omega_m*Omega_p;
        O23[0] = tmp1*(gamma_1[0]+gamma_2[0]-gamma_3[0]);
        O23[1] = tmp1*(gamma_1[1]+gamma_2[1]-gamma_3[1]);

        O32[0] = O23[0];
        O32[1] = O23[1];

        tmp1 = O23[0];
        O23[0] = (eta[0]*eta[0]+eta[1]*eta[1])*tmp1+2*eta[0]*eta[1]*O23[1];
        O23[1] = (eta[0]*eta[0]+eta[1]*eta[1])*O23[1]-2*eta[0]*eta[1]*tmp1;

        tmp1 = O32[0];
        O32[0] = (eta[0]*eta[0]-eta[1]*eta[1])*tmp1-2*eta[0]*eta[1]*O32[1];
        O32[1] = (eta[0]*eta[0]-eta[1]*eta[1])*O23[1]+2*eta[0]*eta[1]*tmp1;

        gamma_1[0] = Psi_1[l][0];
        gamma_1[1] = Psi_1[l][1];
        gamma_2[0] = Psi_2[l][0];
        gamma_2[1] = Psi_2[l][1];
        gamma_3[0] = Psi_3[l][0];
        gamma_3[1] = Psi_3[l][1];

        Psi_1[l][0] = (O11[0]*gamma_1[0]-O11[1]*gamma_1[1]) + (O12[0]*gamma_2[0]-O12[1]*gamma_2[1]) + (O13[0]*gamma_3[0]-O13[1]*gamma_3[1]);
        Psi_1[l][1] = (O11[0]*gamma_1[1]+O11[1]*gamma_1[0]) + (O12[0]*gamma_2[1]+O12[1]*gamma_2[0]) + (O13[0]*gamma_3[1]+O13[1]*gamma_3[0]);

        Psi_2[l][0] = (O21[0]*gamma_1[0]-O21[1]*gamma_1[1]) + (O22[0]*gamma_2[0]-O22[1]*gamma_2[1]) + (O23[0]*gamma_3[0]-O23[1]*gamma_3[1]);
        Psi_2[l][1] = (O21[0]*gamma_1[1]+O21[1]*gamma_1[0]) + (O22[0]*gamma_2[1]+O22[1]*gamma_2[0]) + (O23[0]*gamma_3[1]+O23[1]*gamma_3[0]);

        Psi_3[l][0] = (O31[0]*gamma_1[0]-O31[1]*gamma_1[1]) + (O32[0]*gamma_2[0]-O32[1]*gamma_2[1]) + (O33[0]*gamma_3[0]-O33[1]*gamma_3[1]);
        Psi_3[l][1] = (O31[0]*gamma_1[1]+O31[1]*gamma_1[0]) + (O32[0]*gamma_2[1]+O32[1]*gamma_2[0]) + (O33[0]*gamma_3[1]+O33[1]*gamma_3[0]);
      }

      //MTime.exit_section("Do_NL_Step_light_ana");
    }
  }
}

int main( int argc, char *argv[] )
{
  if ( argc != 2 )
  {
    printf( "No parameter xml file specified.\n" );
    return EXIT_FAILURE;
  }

  ParameterHandler params(argv[1]);
  int dim=0;

  try
  {
    std::string tmp = params.Get_simulation("DIM");
    dim = std::stod(tmp);
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }

  MPI_Init(&argc, &argv);
  fftw_mpi_init();

  try
  {
    if ( dim == 2 )
    {
      MPI::RT_Solver_2::Bragg_double<MPI::Fourier::cft_2d_MPI,2> rtsol( &params );
      rtsol.run_sequence();
    }
    else if ( dim == 3 )
    {
      MPI::RT_Solver_2::Bragg_double<MPI::Fourier::cft_3d_MPI,3> rtsol( &params );
      rtsol.run_sequence();
    }
    else
    {
      cout << "You have found a new dimension!" << endl;
    }
  }
  catch (mu::Parser::exception_type &e)
  {
    cout << "Message:  " << e.GetMsg() << "\n";
    cout << "Formula:  " << e.GetExpr() << "\n";
    cout << "Token:    " << e.GetToken() << "\n";
    cout << "Position: " << e.GetPos() << "\n";
    cout << "Errc:     " << e.GetCode() << "\n";
  }
  catch (std::string &str)
  {
    cout << str << endl;
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
