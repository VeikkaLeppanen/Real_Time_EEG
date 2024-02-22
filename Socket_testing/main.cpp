#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <csignal>
#include <vector>
#include "SamplePacket.h"

volatile std::sig_atomic_t signal_received = 0;

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
        std::cerr << "Socket creation failed" << '\n';
        exit(EXIT_FAILURE);
    }

    memset(&servaddr, 0, sizeof(servaddr));
    memset(&cliaddr, 0, sizeof(cliaddr));
    
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = INADDR_ANY;
    servaddr.sin_port = htons(PORT);
    
    if (bind(sockfd, (const struct sockaddr *)&servaddr, sizeof(servaddr)) < 0) {
        std::cerr << "Bind failed" << '\n';
        close(sockfd); // Ensure the socket is closed on failure
        exit(EXIT_FAILURE);
    }
  
    // Set socket buffer size to 10 MB
    int buffer_size = 1024 * 1024 * 10;
    if (setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, &buffer_size, sizeof(buffer_size)) == -1) {
        std::cerr << "Failed to set socket receive buffer size." << std::endl;
    }

    unsigned char buffer[BUFFER_SIZE];
    socklen_t len = sizeof(cliaddr);
    while (!signal_received) {
        int n = recvfrom(sockfd, (char*)buffer, BUFFER_SIZE, MSG_WAITALL, (struct sockaddr*)&cliaddr, &len);
        if (n < 0) {
            if (signal_received) break; // Check if the signal caused recvfrom to fail
            std::cerr << "Receive failed" << '\n';
            break; // Optionally break on other errors too
        }

        unsigned char firstByte = buffer[0];
        // Handle packets
        switch (firstByte)
        {
        case 0x01: { // MeasurementStart
            /* code */
            break;

        } case 0x02: { // Samples

            // Convert the received data to a vector (if using the vector-based version)
            std::vector<uint8_t> receivedData(buffer, buffer + n);

            // Deserialize the received data into a sample_packet instance
            sample_packet packet = deserializeSamplePacket(receivedData);
            // OR, if using the pointer and size version:
            // sample_packet packet = deserializeSamplePacket(buffer, n);

            // Here you can now access the fields of `packet` directly,
            // for example, print out some deserialized values for verification
            std::cout << "Deserialized PacketSeqNo: " << packet.PacketSeqNo << std::endl;
            std::cout << "NumChannels: " << packet.NumChannels << ", NumSampleBundles: " << packet.NumSampleBundles << std::endl;

            // If you need to process or print the sample data, do it here
            printSamplePacket(packet);

            break;
        } case 0x03: { // Trigger
            /* code */
            break;

        } case 0x04: { // MeasurementEnd
            /* code */
            break;
        
        } case 0x05: { // HardwareState
            /* code */
            break;
        
        
        } default:
            break;
        }
        // Convert the received data to a vector (if using the vector-based version)
        std::vector<uint8_t> receivedData(buffer, buffer + n);

        // Deserialize the received data into a sample_packet instance
        sample_packet packet = deserializeSamplePacket(receivedData);
        // OR, if using the pointer and size version:
        // sample_packet packet = deserializeSamplePacket(buffer, n);

        // Here you can now access the fields of `packet` directly,
        // for example, print out some deserialized values for verification
        std::cout << "Deserialized PacketSeqNo: " << packet.PacketSeqNo << std::endl;
        std::cout << "NumChannels: " << packet.NumChannels << ", NumSampleBundles: " << packet.NumSampleBundles << std::endl;

        // If you need to process or print the sample data, do it here
        printSamplePacket(packet);
    }

    std::cout << "Shutting down..." << '\n';
    close(sockfd);
    return 0;
}
