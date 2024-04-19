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

enum EegBridgeStatus {
  WAITING_MEASUREMENT_START,
  MEASUREMENT_IN_PROGRESS
};

class EegBridge {
public:
    EegBridge() {};

    void bind_socket();
    void spin(dataHandler &handler, volatile std::sig_atomic_t &signal_received);
    void setPort(int port) { PORT = port; }
    
    bool isRunning() { return running; }

private:
    bool running = false;
    EegBridgeStatus eeg_bridge_status;
    
    int numChannels;        // Number of total channels
    int numDataChannels;    // Number of channels storing EEG data
    int numBundles;
    int sampling_rate = 5000;
    int delivery_rate = 5000;
    int lastSequenceNumber = -1;
    
    int PORT = 50000;
    int sockfd;
    struct sockaddr_in servaddr, cliaddr;
    unsigned char buffer[BUFFER_LENGTH];
    Eigen::MatrixXd data_handler_samples;
    socklen_t len;
};

#endif // EEGBRIDGE_H
