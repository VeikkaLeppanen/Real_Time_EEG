#!/bin/bash

# rm -r build/CMakeFiles/real_time_eeg.dir/eegwindow.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/glwidget.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/mainglwidget.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/mainwindow.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/processingworker.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/main.cpp.o
# rm -r build

# export CC=/usr/bin/gcc-10
# export CXX=/usr/bin/g++-10

mkdir build
cd build
# export NIBRARY_ROOT="/home/veikka/Work/EEG/DataStream/Real_Time_EEG/external/nibrary/build-static"
export NIBRARY_ROOT="/home/veikka/Work/EEG/nibrary/build-static"
# export NIBRARY_ROOT="/home/veikka/Work/nibrary-dev/build-static"

# Add debug flags
# export CXXFLAGS="-g -O0"
# export CFLAGS="-g -O0"

cmake .. #-DCMAKE_C_COMPILER=/usr/bin/gcc-10 -DCMAKE_CXX_COMPILER=/usr/bin/g++-10
make

# Set capabilities before running gdb
sudo setcap cap_ipc_lock+ep real_time_eeg

# Run gdb with sudo to ensure it has the right permissions
./real_time_eeg