#include "eeg_bridge.h"

const int PORT = 8080;

volatile std::sig_atomic_t signal_received = 0;
void signal_handler(int signal) {
    signal_received = 1;
}

void eegBridge::bind_socket() {

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

    len = sizeof(cliaddr);
}

void eegBridge::spin() {
    while (!signal_received) {
        int n = recvfrom(sockfd, (char*)buffer, BUFFER_LENGTH, MSG_WAITALL, (struct sockaddr*)&cliaddr, &len);
        if (n < 0) {
            if (signal_received) break; // Check if the signal caused recvfrom to fail
            std::cerr << "Receive failed" << '\n';
            break; // Optionally break on other errors too
        }

        unsigned char firstByte = buffer[0];
        // Handle packets
        switch (firstByte)
        {
        case 0x02: { // Samples

            // Convert the received data to a vector (if using the vector-based version)
            // std::vector<uint8_t> receivedData(buffer, buffer + n);

            // Deserialize the received data into a sample_packet instance
            // sample_packet packet = deserializeSamplePacket_vector(receivedData);
            // OR, if using the pointer and size version:
            sample_packet packet_info;
            std::vector<std::vector<double>> sample_data = deserializeSamplePacket_pointer(buffer, n, packet_info);

            // If you need to process or print the sample data, do it here
            printSamplePacket(packet_info);

            for (int i = 0; i < sample_data.size(); i++) {
                for (int j = 0; j < sample_data[i].size(); ++j) {
                    std::cout << sample_data[i][j] << ' ';
                }
                std::cout  << '\n';
            }

            break;
        } case 0x01: { // MeasurementStart
            /* code */
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
    }
    std::cout << "Shutting down..." << '\n';
    close(sockfd);
}