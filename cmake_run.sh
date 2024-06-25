#!/bin/bash

# rm -r build/CMakeFiles/real_time_eeg.dir/eegwindow.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/glwidget.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/mainglwidget.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/mainwindow.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/processingworker.cpp.o
# rm -r build/CMakeFiles/real_time_eeg.dir/main.cpp.o
# rm -r build

mkdir build
cd build
cmake ..
make
# sudo setcap cap_ipc_lock+ep real_time_eeg
sudo ./real_time_eeg


