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
export NIBRARY_ROOT="/home/veikka/Work/EEG/DataStream/Real_Time_EEG/external/nibrary/build-static"
cmake .. #-DCMAKE_C_COMPILER=/usr/bin/gcc-10 -DCMAKE_CXX_COMPILER=/usr/bin/g++-10
make
sudo setcap cap_ipc_lock+ep real_time_eeg
./real_time_eeg


