#!/bin/bash

# rm -r build/CMakeFiles/real_time_eeg.dir/eegwindow.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/glwidget.cpp.o
rm -r build

mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=~/Qt/6.7.0/gcc_64/lib/cmake -DCMAKE_BUILD_TYPE=Release ..
# cmake ..
make
sudo ./real_time_eeg
