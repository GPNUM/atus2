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
#include <cstring>
#include <ostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <omp.h>
#include "my_structs.h"
#include "noise3_2d.h"
#include "ParameterHandler.h"
#include "fftw3.h"

using namespace std;

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
      cout <<  "Try: " << argv[0] << " noise.bin d0_new.bin [NT]" << endl;
      return EXIT_FAILURE;
    }

  generic_header header;
  ifstream fnoise(argv[1], ifstream::binary );
  if (fnoise.fail()) {
    cout << "File not found: " << argv[1] << endl;
    exit(EXIT_FAILURE);
  }
  fnoise.read( (char*)&header, sizeof(generic_header));
  int64_t NX = header.nDimX;
  int64_t NY = header.nDimY;

  int noise_expansion = 8;
  double dt = header.dx;
  if (argc >= 4) {
    NX = atoi(argv[3])/dt/2;
  }
  double t = atof(argv[3]);
  cout << "t: " << t << endl;
  cout << "NX = " << NX << endl;
  assert (NX <= header.nDimX);

  int blocksize = 1024;
  if ((NX % blocksize) != 0) {
    cout << "Warning: Bad Blocksize: " << NX << " " << blocksize << endl;
    if ((NX % 1000) == 0) {
      blocksize = 1000;
    } else {
      blocksize = 100;
    }
  }
  cout << "Using: " << blocksize << endl;

  if ((NX % blocksize) != 0) {
    cout << "Still Bad Blocksize!" << endl;
  }

  int64_t blocks = NX/blocksize;
  cout << "Blocks " << blocks << endl;

  double *noise;
  noise = new double[ blocksize * NY ];
  double *integrated_noise = new double[NY];
  for (int64_t i = 0; i < NY; i++) {
    integrated_noise[i] = 0.0;
  }

  ofstream ftotal( "int_total.txt");
  for (int64_t block = 0; block < blocks; block++) {
    fnoise.read( (char*)noise, sizeof(double)*blocksize*NY );
    #pragma omp parallel for
    for (int64_t j = 0; j < blocksize; j++) {
      for (int64_t i = 0; i < NY; i++) {
        integrated_noise[i] += noise[i+j*NY];
      }
    }

    double total = 0.0;
    for (int64_t i = 0; i < NY; i++) {
      total += integrated_noise[i]*header.dy*header.dx;
    }
    total *= header.dx;
    ftotal << block << "\t" << total << "\t" << integrated_noise[NY/2]*header.dx << endl;

  }

  for (int64_t i = 0; i < NY; i++) {
    integrated_noise[i] *= header.dx*0.5;
  }

  double total = 0.0;
  for (int64_t i = 0; i < NY; i++) {
    total += integrated_noise[i];
  }
  total *= header.dy;
  cout << "Total\t" << total << endl;

  ifstream fd0_new(argv[2], ifstream::binary );
  if (fd0_new.fail()) {
    cout << "File not found: " << argv[2] << endl;
    exit(EXIT_FAILURE);
  }
  fd0_new.seekg(sizeof(generic_header));
  double *d0_new;
  d0_new = new double[ blocksize * NY ];
  double *integrated_d0_new = new double[NY];
  for (int64_t i = 0; i < NY; i++) {
    integrated_d0_new[i] = 0.0;
  }

  ofstream ftotal_d0_new( "int_total_d0_new.txt");
  for (int64_t block = 0; block < blocks; block++) {
    fd0_new.read( (char*)d0_new, sizeof(double)*blocksize*NY );
    #pragma omp parallel for
    for (int64_t j = 0; j < blocksize; j++) {
      for (int64_t i = 0; i < NY; i++) {
        integrated_d0_new[i] += d0_new[i+j*NY];
      }
    }

    double total_d0_new = 0.0;
    for (int64_t i = 0; i < NY; i++) {
      total_d0_new += integrated_d0_new[i];
    }
    total_d0_new *= header.dx;
    ftotal_d0_new << block << "\t" << total_d0_new << endl;

  }

  for (int64_t i = 0; i < NY; i++) {
    integrated_d0_new[i] *= header.dx*0.5;
  }

  double total_d0_new = 0.0;
  for (int64_t i = 0; i < NY; i++) {
    total_d0_new += integrated_d0_new[i];
  }
  total_d0_new *= header.dy;
  cout << "Total_D0_New\t" << total_d0_new << endl;

  ifstream fd2(argv[4], ifstream::binary);
  ifstream fd1(argv[5], ifstream::binary);
  if (fd2.fail()) {
    cout << "File not found: " << argv[4] << endl;
    exit(EXIT_FAILURE);
  }
  if (fd1.fail()) {
    cout << "File not found: " << argv[5] << endl;
    exit(EXIT_FAILURE);
  }
  double *d2 = new double[NY];
  double *d1 = new double[NY];
  fd2.seekg(sizeof(generic_header) + sizeof(double)*(NX-1)*NY);
  fd1.seekg(sizeof(generic_header) + sizeof(double)*(NX-1)*NY);
  fd2.read( (char*)d2, sizeof(double)*NY );
  fd1.read( (char*)d1, sizeof(double)*NY );

  //  double variance = 427.68;
  double variance = 0.5;
  double sigma = sqrt(variance);
  double noise_strength = 0.5;
  double epsilon = noise_strength;
  double g = 1.0/(2.0*variance);
  double *gauss = new double[NY];
  double N = 10000.0;
  double garg = 0.0;
  for (int64_t i = 0; i < NY; ++i) {
    garg += -d1[i]/(2*(1-d2[i]));
  }
  cout << "Garg: " << garg << " G: " << exp(garg) << endl;
  cout << "Sigma: " << sigma << endl;


  double alpha = 0.000365;
  double alpha2 = pow(alpha,2);
  for (int64_t i = 0; i < NY; i++) {
    double A = integrated_noise[i];
    A  =0;
    double x = header.yMin+i*header.dy;
    double gamma = sqrt((pow(sigma,2) + sqrt(pow(sigma,4)+pow(A,2)))/2);
    garg = 0.0;
    double phase = integrated_d0_new[i];
    //    A = 0;
    //    phase = 0;
    //gauss[i] = N/sqrt(2.0*M_PI)/sigma*exp(-pow(x,2)/(2*pow(sigma,2))); // Without Noise
    //gauss[i] = N/sqrt(2.0*M_PI)*sigma/(sqrt(pow(sigma,4)+pow(A,2)))*exp(-pow(x,2)*pow(sigma,2)/(2*(pow(sigma,4) + pow(A,2)))); // 0. Order
    //    gauss[i] *= exp(2*garg);

    //    gauss[i] += pow(epsilon,2)*N/sqrt(2.0*M_PI)/sigma*exp(-pow(x,2)/(2.0*pow(sigma,2)));
     // gauss[i] += 0.01*epsilon*2*N/sqrt(2.0*M_PI)/(sqrt(pow(sigma,4)+pow(A,2)))*exp(-pow(x,2)/4.0*((pow(sigma,2)/(pow(sigma,4)+pow(A,2)) + 1.0/pow(sigma,2))))*(sqrt(0.5*(pow(sigma,2) + sqrt(pow(sigma,4) + pow(A,2)))) * cos(pow(x,2)/4*A/(pow(sigma,4)+ pow(A,2)) + phase)
     //                                                                                                                                               + sqrt(0.5*(-pow(sigma,2) + sqrt(pow(sigma,4) + pow(A,2)))) * sin(pow(x,2)/4*A/(pow(sigma,4)+ pow(A,2)) + phase));

    //    gauss[i] = N/sqrt(2.0*M_PI)*sigma/sqrt(pow(sigma,4) + alpha2*pow(t,2))*exp(-pow(x,2)*pow(sigma,2)/(2*(pow(sigma,4) + alpha2*pow(t,2)))); // With time Without Noise

    // With Time
    gauss[i] = N/sqrt(2.0*M_PI)*sigma/(sqrt(pow(sigma,4)+alpha2*pow(t+A,2)))*exp(-pow(x,2)*pow(sigma,2)/(2*(pow(sigma,4) + alpha2*pow(t+A,2)))); // 0. Order

    double a = pow(sigma,4) + alpha2 * t*(t+A);
    double b = pow(sigma,2)*alpha*A;
    //    gauss[i] += pow(epsilon,2)*N/sqrt(2.0*M_PI)*sigma/sqrt(pow(sigma,4) + alpha2*pow(t,2))*exp(-pow(x,2)*pow(sigma,2)/(2.0*(pow(sigma,4) + alpha2*pow(t,2))));
    //    gauss[i] += epsilon*2*N/sqrt(2.0*M_PI)*sigma/(sqrt(pow(sigma,4)+alpha2*pow(t+A,2))*sqrt(pow(sigma,4)+alpha2*pow(t,2)))*exp(-pow(x,2)/4*(pow(sigma,2)/(pow(sigma,4) + alpha2*pow(t+A,2)) + pow(sigma,2)/(pow(sigma,4) + alpha2*pow(t,2))))
    //      * (sqrt(0.5*(a+sqrt(pow(a,2)+pow(b,2))))*cos(pow(x,2)*alpha/4*(t/(pow(sigma,4) + alpha2*pow(t,2)) - (t+A)/(pow(sigma,4) + alpha2*pow(t+A,2))) - alpha*phase)
    //         + sqrt(0.5*(-a+sqrt(pow(a,2)+pow(b,2))))*sin(pow(x,2)*alpha/4*(t/(pow(sigma,4) + alpha2*pow(t,2)) - (t+A)/(pow(sigma,4) + alpha2*pow(t+A,2))) - alpha*phase));

    //    gauss[i] *= exp(2*garg);
 }

  header.nDims = 1;
  header.nDimX = header.nDimY;
  header.nDimY = 1;
  header.nDimZ = 1;
  header.dx = header.dy;
  header.xMin = header.yMin;
  header.xMax = header.yMax;

  for (int64_t i = 0; i < NY; ++i) {
    integrated_noise[i] *= integrated_noise[i];
  }

  char* bin_header = reinterpret_cast<char*>(&header);
  string foo = string(argv[1]);
  ofstream file1( "int_"+foo, ofstream::binary );
  file1.write( bin_header, sizeof(generic_header) );
  char* bin_noise;
  bin_noise = reinterpret_cast<char*>(integrated_noise);
  file1.write( bin_noise, NY*sizeof(double) );

  ofstream file2( "noise_gauss.bin", ofstream::binary );
  file2.write( bin_header, sizeof(generic_header) );
  char* bin_gauss;
  bin_gauss = reinterpret_cast<char*>(gauss);
  file2.write( bin_gauss, NY*sizeof(double) );
  file2.close();


  for (int64_t i = 0; i < NY; ++i) {
    double A = 0.0;
    double x = header.yMin+i*header.dy;
    gauss[i] = N/sqrt(2.0*M_PI)*sigma/(sqrt(pow(sigma,4)+alpha2*pow(t+A,2)))*exp(-pow(x,2)*pow(sigma,2)/(2*(pow(sigma,4) + alpha2*pow(t+A,2))));
  }
  ofstream file3("gauss.bin", ofstream::binary);
  file3.write(bin_header, sizeof(generic_header) );
  file3.write( bin_gauss, NY*sizeof(double) );
  file3.close();

  file1.close();
  fnoise.close();
  delete[] noise;
  delete[] integrated_noise;

}
