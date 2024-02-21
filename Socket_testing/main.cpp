#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <csignal>
#include <vector>

volatile std::sig_atomic_t signal_received = 0;

struct sample_packet {
    uint8_t FrameType;
    uint8_t MainUnitNum;
    uint8_t Reserved[2];
    uint32_t PacketSeqNo;
    uint16_t NumChannels;
    uint16_t NumSampleBundles;
    uint64_t FirstSampleIndex;
    uint64_t FirstSampleTime;
    uint32_t Samples[1][4];
};

std::vector<std::vector<uint32_t>> process24BitSamples(const std::vector<uint8_t>& rawData, uint16_t numChannels, uint16_t numSampleBundles) {
    std::vector<std::vector<uint32_t>> samples(numChannels, std::vector<uint32_t>(numSampleBundles));

    size_t index = 0; // This would start after the sample_packet data in rawData
    for (uint16_t bundle = 0; bundle < numSampleBundles; ++bundle) {
        for (uint16_t channel = 0; channel < numChannels; ++channel) {
            if (index + 3 <= rawData.size()) {
                // Assuming samples are stored in little-endian format
                samples[channel][bundle] = (rawData[index] | (rawData[index + 1] << 8) | (rawData[index + 2] << 16));
                index += 3;
            }
        }
    }

    return samples;
}


// Temporary termination method for the receiving loop with Ctrl+C
void signal_handler(int signal) {
    signal_received = 1;
}

const int PORT = 8080;
const int BUFFER_SIZE = 1472;

int main() {
    int sockfd;
    struct sockaddr_in servaddr, cliaddr;

    // Register signal handler
    std::signal(SIGINT, signal_handler);

    sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0) {
        std::cerr << "Socket creation failed" << std::endl;
        exit(EXIT_FAILURE);
    }

    memset(&servaddr, 0, sizeof(servaddr));
    memset(&cliaddr, 0, sizeof(cliaddr));
    
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = INADDR_ANY;
    servaddr.sin_port = htons(PORT);
    
    if (bind(sockfd, (const struct sockaddr *)&servaddr, sizeof(servaddr)) < 0) {
        std::cerr << "Bind failed" << std::endl;
        close(sockfd); // Ensure the socket is closed on failure
        exit(EXIT_FAILURE);
    }

    unsigned char buffer[BUFFER_SIZE];
    socklen_t len = sizeof(cliaddr);

    while (!signal_received) {
        int n = recvfrom(sockfd, (char *)buffer, BUFFER_SIZE, MSG_WAITALL, (struct sockaddr *) &cliaddr, &len);
        if (n < 0) {
            if (signal_received) break; // Check if the signal caused recvfrom to fail
            std::cerr << "Receive failed" << std::endl;
            break; // Optionally break on other errors too
        }
        buffer[n] = '\0';
        std::cout << "Received: " << std::size(buffer) << std::endl;
    }

    std::cout << "Shutting down..." << std::endl;
    close(sockfd);
    return 0;
}
