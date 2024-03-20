#!/bin/bash
g++ -fopenmp $(pkg-config --cflags python3) -o build/datastream main.cpp eeg_bridge/eeg_bridge.cpp eeg_bridge/measurementStartPacket.cpp eeg_bridge/samplePacket.cpp eeg_bridge/networkUtils.cpp dataHandler/dataHandler.cpp dataHandler/GACorrection.cpp -lpython3.11 -I/usr/include/eigen3 -I/usr/lib/python3.11/site-packages/numpy/core/include -Wno-deprecated-declarations

./build/datastream