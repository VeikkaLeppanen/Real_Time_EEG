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
#include <boost/stacktrace.hpp>

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
    void spin(volatile std::sig_atomic_t &signal_received);
    void setPort(int port) { PORT = port; }
    void setTimeout(int timeout) { socket_timeout = timeout; }
    
    int receive_packet() { return recvfrom(sockfd, (char*)buffer, BUFFER_LENGTH, MSG_WAITALL, (struct sockaddr*)&cliaddr, &len); }
    void close_socket() { close(sockfd); }

    bool isRunning() { return running; }

    bool running = false;
    EegBridgeStatus eeg_bridge_status;
    
    int numChannels;        // Number of total channels
    int numDataChannels;    // Number of channels storing EEG data
    int numBundles;
    int sampling_rate = 5000;
    int delivery_rate = 5000;
    int lastSequenceNumber = -1;

    const uint8_t DC_MODE_SCALE = 100;
    const uint16_t NANO_TO_MICRO_CONVERSION = 1000;
    const double DOUBLESCALINGFACTOR = 10000.0;

    unsigned char buffer[BUFFER_LENGTH];
    Eigen::MatrixXd data_handler_samples;

private:
    int PORT = 50000;
    int socket_timeout = 60;
    int sockfd;
    struct sockaddr_in servaddr, cliaddr;
    socklen_t len;
};

#endif // EEGBRIDGE_H
