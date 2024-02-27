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

// The maximum length of the UDP packet, as mentioned in the manual of Bittium NeurOne.
#define BUFFER_LENGTH 1472

class eegBridge {
public:
    eegBridge() {};

    void bind_socket();
    void spin();

private:
    int numChannels = 4;
    int numBundles = 2;
    int sockfd;
    struct sockaddr_in servaddr, cliaddr;
    unsigned char buffer[BUFFER_LENGTH];
    socklen_t len;
};

#endif // EEGBRIDGE_H
