#ifndef EEGBRIDGE_H
#define EEGBRIDGE_H

#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <csignal>
#include <vector>
#include "samplePacket.h"
#include "measurementStartPacket.h"
#include "../dataHandler/dataHandler.h"

// The maximum length of the UDP packet, as mentioned in the manual of Bittium NeurOne.
#define BUFFER_LENGTH 1472

enum EegBridgeState {
  WAITING_FOR_MEASUREMENT_START,
  WAITING_FOR_MEASUREMENT_STOP,
  WAITING_FOR_SESSION_STOP,
  WAITING_FOR_SESSION_START,
  STREAMING,
  ERROR_OUT_OF_SYNC,
  ERROR_SAMPLES_DROPPED
};

class EegBridge {
public:
    EegBridge() {};

    void bind_socket();
    void spin(dataHandler &handler);

private:
    EegBridgeState eeg_bridge_state;
    
    int numChannels;        // Number of total channels
    int numDataChannels;    // Number of channels storing EEG data
    int numBundles;
    int sampling_rate = 5000;
    int delivery_rate = 5000;
    int packet_sequence_number = -1;
    
    int sockfd;
    struct sockaddr_in servaddr, cliaddr;
    unsigned char buffer[BUFFER_LENGTH];
    Eigen::MatrixXd data_handler_samples;
    socklen_t len;
};

#endif // EEGBRIDGE_H
