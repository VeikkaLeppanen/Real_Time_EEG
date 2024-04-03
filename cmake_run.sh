#!/bin/bash

rm -r build
mkdir build
cd build
cmake ..
make
./real_time_eeg
