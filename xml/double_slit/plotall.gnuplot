#!/usr/bin/gnuplot -persist
#
set pm3d map

LIST_PNG=system("export LANG=us_US.UTF-8 && seq -f '%.3f_1.jpg'  0.1 0.1 2")
LIST_BIN=system("export LANG=us_US.UTF-8 && seq -f '\"<gpo3 %.3f_1.bin\"'  0.1 0.1 2")

set terminal jpeg enhanced size 1024,1024

do for [i=1:20:1] {
  set output word(LIST_PNG,i)
  splot word(LIST_BIN,i)
}
#    EOF
