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

#include "timer.h"

void Timer::write2file()
{
  FILE *fh = nullptr;

  fh = fopen( "Process_time.txt", "w" );
  for ( map<string,elements>::iterator it = function.begin(); it != function.end(); it++)
  {
    fprintf(fh,"%s:\nCounter: %d\nAbsCounter: %d\nWalltime: %g\nRealtime: %g\n",it->first.c_str(),it->second.counter,it->second.abs_counter,it->second.walltime,it->second.realtime);
    fprintf(fh,"--------------------------------------------------------------------------------------\n");
  }
  fclose(fh);
}

void Timer::write2file2()
{
  FILE *fh = nullptr;

  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( rank == 0)
  {
    if ( m_counter == 1) fh = fopen( "Process_time.txt", "w" );
    else fh = fopen( "Process_time.txt", "a" );

    for ( map<string,elements>::iterator it = function.begin(); it != function.end(); it++)
    {
      fprintf(fh,"Durchlauf %d:\n",m_counter);
      fprintf(fh,"%s:\nCounter: %d\nAbsCounter: %d\nRealtime: %g\nCPU-Time: %g\n",it->first.c_str(),it->second.counter,it->second.abs_counter,it->second.walltime,it->second.realtime);
      fprintf(fh,"--------------------------------------------------------------------------------------\n");
    }

    fclose(fh);
  }
}

void Timer::enter_section(string s)
{
  function[s].starttime = MPI_Wtime();
  function[s].counter++;
}

void Timer::exit_section(string s)
{
  function[s].endtime = MPI_Wtime();
  function[s].walltime += function[s].endtime-function[s].starttime;
  MPI_Allreduce(&function[s].walltime,&function[s].realtime,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&function[s].counter,&function[s].abs_counter,1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}
