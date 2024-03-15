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

class EegBridge {
public:
    EegBridge() {};

    void bind_socket();
    void spin(dataHandler &handler);

private:
    int numChannels = 4;
    int numBundles = 2;
    int sampling_rate = 5000;
    int delivery_rate = 5000;
    int sockfd;
    struct sockaddr_in servaddr, cliaddr;
    unsigned char buffer[BUFFER_LENGTH];
    Eigen::MatrixXd data_handler_samples;
    socklen_t len;
};

#endif // EEGBRIDGE_H
