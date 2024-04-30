#!/bin/bash

# rm -r build/CMakeFiles/real_time_eeg.dir/main.cpp.o
rm -r build

mkdir build
cd build
cmake ..
make
./real_time_eeg
