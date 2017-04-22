#!/bin/sh
mkdir -p $HOME/local/src/
cd $HOME/local/src/
git clone -b release https://github.com/roswell/roswell.git
cd roswell
sh bootstrap
./configure --prefix=$HOME/local/
make
make install
ros setup
