#!/bin/sh
mkdir -p $HOME/.local/share/common-lisp/systems/
ln -s $(pwd) $HOME/.local/share/common-lisp/systems/atus2

mkdir -p $HOME/.roswell/local-projects/$USER/
ln -s $(pwd) $HOME/.roswell/local-projects/$USER/atus2
