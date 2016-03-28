#!/bin/bash
ctags -R --fields=+lS --c++-kinds=+p --fields=+iaS --extra=+q .
g++ -std=c++0x -o myprogram.exe -static -static-libgcc -static-libstdc++ *.cpp

