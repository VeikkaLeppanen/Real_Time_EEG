#!/bin/bash

# rm -r build/CMakeFiles/real_time_eeg.dir/eegwindow.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/glwidget.cpp.o
# rm -r build

mkdir build
cd build
export NIBRARY_ROOT="/home/user/nibrary/build-static"
cmake -DCMAKE_PREFIX_PATH=~/Qt/6.7.0/gcc_64/lib/cmake -DCMAKE_BUILD_TYPE=Release ..
# cmake ..
make # VERBOSE=1
sudo setcap cap_ipc_lock+ep real_time_eeg
sudo sudo chmod 666 /dev/ttyUSB0
./real_time_eeg
