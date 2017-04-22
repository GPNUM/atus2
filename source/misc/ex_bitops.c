/* * ATUS2 - The ATUS2 package is atom interferometer Toolbox developed at ZARM
 * (CENTER OF APPLIED SPACE TECHNOLOGY AND MICROGRAVITY), Germany. This project is
 * founded by the DLR Agentur (Deutsche Luft und Raumfahrt Agentur). Grant numbers:
 * 50WM0942, 50WM1042, 50WM1342.
 * Copyright (C) 2017 Želimir Marojević, Ertan Göklü, Claus Lämmerzahl
 *
 * This file is part of ATUS2.
 *
 * ATUS2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ATUS2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ATUS2.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
./arguments -a -b
./arguments -abc
./arguments -h
./arguments -z
*/

#include <stdio.h>

//anonymous enum of options, can be done any way you wish
enum {
    OPT_A = 0x01,
    OPT_B = 0x02,
    OPT_C = 0x04,
    OPT_H = 0x08,
};

int main(int argc, char **argv)
{
    //unsigned int, or rather uint32_t = 8*4=32 bits for options if needed, unsigned means last bit is ours to use
    unsigned int opt = 0x0;
    //can do char array for options like '-nodebug'
    char c;

    //extract arguments from argument array.
    while((++argv)[0] && argv[0][0] == '-')
    {
        while((c = *++argv[0]) != 0)
        {
            switch(c) { 
            case 'a':
            //assign option bits to "opt" bit array
                opt |= OPT_A;  break;
            case 'b':
                opt |= OPT_B;  break;
            case 'c':
                opt |= OPT_C;  break;
            case 'h':
                opt |= OPT_H;  break;
            //this will happen if they enter an invalid option:
            default: 
                printf("%s: Unknown option %c", argv[0], c);
                return 1; //break out of application
            }
        }
    }
    
    //apply bitwise AND to check for assignedness a few times
    if(opt & OPT_A)
        printf("Hello World!\n");
   
    if(opt & OPT_B) {
        unsigned int foo;
        foo = 2000;
        printf("Foo has been initialized.\n");
    }
    
    //compare if two flags were specifically set
    if ((opt & (OPT_B | OPT_C)) == (OPT_B | OPT_C)) 
        printf("Flags B and C were set.\n");
    
    if(opt & OPT_H) {
        //print help, may wish to create exit point to stop program from executing
        printf("\tHelp is not implemented yet\n\tAllowable options: [-abch]\n");
        return 0;
    }
  
    //----------------- Some fun extras: ---------------------//
    
    //Reset bitflag completely
    opt = 0;
    
    //Apply bitwise OR to append multiple flags
    opt = (OPT_A | OPT_B | OPT_C);
    
    //Apply bitwise AND+EQUALS to add or remove flags to existing option field
    //Then we apply bitwise NOT (a complement) to remove both flags
    opt &= ~(OPT_A | OPT_B);
    
    //Options A and B are now removed
    
    //Check if BOTH flags are not set
    if ((opt & (OPT_A | OPT_B)) == 0)
        //printf( A and B are not set )

    //check if only one is not set
    if ((opt & OPT_A) == 0)
        //printf( flag A is not set )
        
    //end program
    return 0;
    
}