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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

int detuning[] = {1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000};

int main( int argc, char *argv[] )
{
  string base_folder = "Rb_K";
  string folder;
  string shellcmd;
  string filename;
  string tmp;

  shellcmd = "mkdir " + base_folder;
  system( shellcmd.c_str() );
  
  vector<string> elem_strings;
  vector<string> text_strings;
   
  filename = base_folder + "/run_all";
  ofstream bash_script(filename);
  bash_script << "#!/bin/bash" << "\n";

  int dN1=10, dN2=10;
  int N1, N2;

  for( int i1=1; i1<50; i1++ )
  {  
    N1 = i1*dN1;
    N2 = i1*dN2;

    for( int j1=0; j1<sizeof(detuning)/sizeof(int); j1++ )
    {
      folder = base_folder + "/" + to_string(N1) + "_" + to_string(N2) + "/" + to_string(detuning[j1]);
      shellcmd = "mkdir -p " + folder; 
      system( shellcmd.c_str() );

      filename = folder + "/params.xml";

      FILE * fh = fopen( filename.c_str(), "w" );

      XMLPrinter printer(fh);
      printer.OpenElement( "SIMULATION" );

        printer.OpenElement( "DIM" );
          printer.PushText("1");
        printer.CloseElement();

        printer.OpenElement( "FILENAME" );
          printer.PushText("../inf_Rb_0.000_1.bin");
        printer.CloseElement();

        printer.OpenElement( "FILENAME_2" );
          printer.PushText("../inf_zero.bin");
        printer.CloseElement();

        printer.OpenElement( "FILENAME_3" );
          printer.PushText("../inf_K_0.000_1.bin");
        printer.CloseElement();

        printer.OpenElement( "FILENAME_4" );
          printer.PushText("../inf_zero.bin");
        printer.CloseElement();

        elem_strings.clear();
        text_strings.clear();
        elem_strings.push_back("Beta"); text_strings.push_back("0");
        elem_strings.push_back("laser_k"); text_strings.push_back("8.05289");
        elem_strings.push_back("laser_k_2"); text_strings.push_back("8.1951");
        elem_strings.push_back("laser_domh"); text_strings.push_back("0.0471239");
        elem_strings.push_back("laser_domh_2"); text_strings.push_back("0.109465");
        elem_strings.push_back("laser_dk"); text_strings.push_back("0");
        elem_strings.push_back("laser_dk_2"); text_strings.push_back("0");
        elem_strings.push_back("Phase"); text_strings.push_back("0");
        elem_strings.push_back("Phase_2"); text_strings.push_back("0");
        elem_strings.push_back("sweep"); text_strings.push_back("0");
        elem_strings.push_back("rabi_threshold"); text_strings.push_back("7");
        elem_strings.push_back("stepsize"); text_strings.push_back("0.001");

        // problem related parameter
        printer.OpenElement( "CONSTANTS" );
          for( int i=0; i<elem_strings.size(); i++ )
          {
            printer.OpenElement( elem_strings[i].c_str() );
              printer.PushText( text_strings[i].c_str() );
            printer.CloseElement();
          }
        printer.CloseElement(); // close CONSTANTS
          
        elem_strings.clear();
        text_strings.clear();
        elem_strings.push_back("Amp_1"); text_strings.push_back("-5.60499,-5.60499");
        elem_strings.push_back("Amp_2"); text_strings.push_back("-10.182,-10.182");
        elem_strings.push_back("Alpha_1"); text_strings.push_back("0.000365368,0.000365368,0.000365368");
        elem_strings.push_back("Alpha_2"); text_strings.push_back("0.00081496,0.00081496,0.00081496");
        //elem_strings.push_back("Delta_L"); text_strings.push_back("0,-6283.19,0,-3300");
        string tmp_str = "0,-" + to_string(detuning[j1]) + ",0,-3300";
        elem_strings.push_back("Delta_L"); text_strings.push_back(tmp_str.c_str());
        //elem_strings.push_back("GS_1"); text_strings.push_back("4.38015e-5,0,2.43e-5,0");
        elem_strings.push_back("GS_1"); text_strings.push_back("0,0,0,0");
        elem_strings.push_back("GS_2"); text_strings.push_back("0,0,0,0");
        //elem_strings.push_back("GS_3"); text_strings.push_back("2.43e-5,0,1.51e-4,0");
        elem_strings.push_back("GS_3"); text_strings.push_back("0,0,0,0");
        elem_strings.push_back("GS_4"); text_strings.push_back("0,0,0,0");
        elem_strings.push_back("Beta"); text_strings.push_back("0,0,0");
        elem_strings.push_back("Beta_2"); text_strings.push_back("0,0,0");
        elem_strings.push_back("P0"); text_strings.push_back("0,0,0");
        string N_str = to_string(N1) + "," + to_string(N2);
        elem_strings.push_back("N"); text_strings.push_back(N_str);
      
        // mesh related parameter
        printer.OpenElement( "VCONSTANTS" );
          for( int i=0; i<elem_strings.size(); i++ )
          {
            printer.OpenElement( elem_strings[i].c_str() );
              printer.PushText( text_strings[i].c_str() );
            printer.CloseElement();
          }
        printer.CloseElement(); // close VCONSTANTS

        // algorithm related parameter
        printer.OpenElement( "AI_SEQUENCE" );
          printer.OpenElement( "bragg_ad" );
            printer.PushAttribute("dt","0.1");
            printer.PushAttribute("Nk","10");
            printer.PushAttribute("output_freq","last");
            printer.PushAttribute("pn_freq","each");
            printer.PushAttribute("rabi_output_freq","each");
            printer.PushAttribute("custom_freq","none");
            printer.PushText( "2000,2000" );
          printer.CloseElement();
        printer.CloseElement(); // close AI_SEQUENCE
      printer.CloseElement(); // close SIMULATION

      if( j1 == 0 )
      {
        filename = base_folder + "/" + to_string(N1) + "_" + to_string(N2) + "/breed.xml";  
        FILE * fh2 = fopen( filename.c_str(), "w" );

        XMLPrinter printer2(fh2);
        printer2.OpenElement( "SIMULATION" );
          
          printer2.OpenElement( "DIM" );
            printer2.PushText("1");
          printer2.CloseElement();

          printer2.OpenElement( "DTYP" );
            printer2.PushText("COMPLEX");
          printer2.CloseElement();

          printer2.OpenElement( "FILENAME" );
            printer2.PushText("Rb_0.000_1.bin");
          printer2.CloseElement();

          printer2.OpenElement( "FILENAME_3" );
            printer2.PushText("K_0.000_1.bin");
          printer2.CloseElement();

          printer2.OpenElement( "GUESS_1D" );
            printer2.PushText("exp(-0.5*x^2)");
          printer2.CloseElement();

          printer2.OpenElement( "GUESS_1D_2" );
            printer2.PushText("x*exp(-0.25*x^2)");
          printer2.CloseElement();

          printer2.OpenElement( "POTENTIAL" );
            //printer2.PushText("1.8727*x^2"); // 1 kHz Falle 
            printer2.PushText("7.4909*x^2"); // 2 kHz Falle 
          printer2.CloseElement();

          printer2.OpenElement( "POTENTIAL_2" );
            //printer2.PushText("0.3764*x^2");  // 1 kHz Falle 
            printer2.PushText("1.5056*x^2");  // 2 kHz Falle 
          printer2.CloseElement();

          // problem related parameter
            
          elem_strings.clear();
          text_strings.clear();
          elem_strings.push_back("Alpha_1"); text_strings.push_back("1,1,1,1");
          elem_strings.push_back("Alpha_2"); text_strings.push_back("1,1,1,1");
          elem_strings.push_back("GS_1"); text_strings.push_back("0.1197,0.0666");
          elem_strings.push_back("GS_2"); text_strings.push_back("0.0298,0.1861");
          elem_strings.push_back("N"); text_strings.push_back(N_str);
        
          // mesh related parameter
          printer2.OpenElement( "VCONSTANTS" );
            for( int i=0; i<elem_strings.size(); i++ )
            {
              printer2.OpenElement( elem_strings[i].c_str() );
                printer2.PushText( text_strings[i].c_str() );
              printer2.CloseElement();
            }
          printer2.CloseElement(); // close VCONSTANTS
            
          printer2.OpenElement( "CONSTANTS" );
            printer2.OpenElement( "rabi_threshold" );
              printer2.PushText( "7" );
            printer2.CloseElement();
            printer2.OpenElement( "stepsize" );
              printer2.PushText( "0.0001" );
            printer2.CloseElement();
          printer2.CloseElement(); // close CONSTANTS

          elem_strings.clear();
          text_strings.clear();
          //elem_strings.push_back("NX"); text_strings.push_back("1024"); // 1kHz
          elem_strings.push_back("NX"); text_strings.push_back("512"); // 2kHz
          elem_strings.push_back("NY"); text_strings.push_back("256"); 
          elem_strings.push_back("NZ"); text_strings.push_back("256"); 
          //elem_strings.push_back("XMIN"); text_strings.push_back("-20"); // 1kHz
          //elem_strings.push_back("XMAX"); text_strings.push_back("20"); // 1kHz
          elem_strings.push_back("XMIN"); text_strings.push_back("-10"); // 2kHz
          elem_strings.push_back("XMAX"); text_strings.push_back("10"); // 2kHz
          elem_strings.push_back("YMIN"); text_strings.push_back("-40");
          elem_strings.push_back("YMAX"); text_strings.push_back("40");
          elem_strings.push_back("ZMIN"); text_strings.push_back("-40");
          elem_strings.push_back("ZMAX"); text_strings.push_back("40");

          // mesh related parameter
          printer2.OpenElement( "ALGORITHM" );
            for( int i=0; i<elem_strings.size(); i++ )
            {
              printer2.OpenElement( elem_strings[i].c_str() );
                printer2.PushText( text_strings[i].c_str() );
              printer2.CloseElement();
            }
          printer2.CloseElement(); // close ALGORITHM
        printer2.CloseElement(); // close SIMULATION 
      
        tmp = to_string(N1) + "_" + to_string(N2);
        bash_script << "cd " + tmp + "/\n";
        bash_script << "sobmin_1d_2 breed.xml\n";
        //bash_script << "inflate_domain --fak 8 K_0.000_1.bin\n"; // 1kHz
        //bash_script << "inflate_domain --fak 8 Rb_0.000_1.bin\n"; // 1kHz
        //bash_script << "inflate_domain --fak 8 zero.bin\n"; // 1kHz
        bash_script << "inflate_domain --fak 16 K_0.000_1.bin\n"; // 2kHz
        bash_script << "inflate_domain --fak 16 Rb_0.000_1.bin\n"; // 2kHz
        bash_script << "inflate_domain --fak 16 zero.bin\n"; // 2kHz
      }

      bash_script << "cd " << to_string(detuning[j1]) << "\n";
      bash_script << "bragg_ds_1d params.xml\n";
      bash_script << "cd ..\n";

    }
    bash_script << "cd ..\n";
  }

  shellcmd = "chmod +x " + base_folder + "/run_all";
  system(shellcmd.c_str());
return 0;
}
