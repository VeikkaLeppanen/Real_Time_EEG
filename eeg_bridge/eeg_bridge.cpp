#include "eeg_bridge.h"

const int PORT = 8080;

// This is used to terminate the program with Ctrl+C
volatile std::sig_atomic_t signal_received = 0;
void signal_handler(int signal) {
    signal_received = 1;
}

void EegBridge::bind_socket() {

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

void EegBridge::spin(dataHandler &handler) {
    std::cout << "Waiting for packets..." << '\n';
    eeg_bridge_state = WAITING_FOR_MEASUREMENT_START;
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
        case 0x02: { // SamplesPacket

            // Deserialize the received data into a sample_packet instance
            sample_packet packet_info;
            deserializeSamplePacketEigen_pointer(buffer, n, packet_info, data_handler_samples);

            // AMPLIFIER *100
            // NANO TO MICRO /1000

            for (int i = 0; i < packet_info.NumSampleBundles; i++) {
                handler.addData(data_handler_samples.col(i), static_cast<double>(packet_info.FirstSampleTime), 0);
            }

            // std::cout << handler.getDataInOrder(1) << '\n';
            // std::cout << handler.get_buffer_capacity() << '\n';
            break;

        } case 0x01: { // MeasurementStartPacket
            
            std::cout << "MeasurementStart package received!\n";
            measurement_start_packet packet_info;
            std::vector<uint16_t> SourceChannels;
            std::vector<uint8_t> ChannelTypes;
            
            deserializeMeasurementStartPacket_pointer(buffer, n, packet_info, SourceChannels, ChannelTypes);

            numChannels = packet_info.NumChannels;
            sampling_rate = packet_info.SamplingRateHz;

            // TODO: Initialize data_handler_samples
            data_handler_samples = Eigen::MatrixXd::Zero(numChannels, 1000);


            std::cout << "MeasurementStart package processed!\n";

            handler.reset_handler(numChannels, sampling_rate);
            std::cout << "DataHandler reset!\n";
            eeg_bridge_state = WAITING_FOR_MEASUREMENT_STOP;

            break;

        } case 0x03: { // TriggerPacket
            /* code */
            break;

        } case 0x04: { // MeasurementEndPacket
            /* code */
            break;
        
        } case 0x05: { // HardwareStatePacket
            /* code */
            break;
        
        
        } default:
            break;
        }
    }
    std::cout << "Shutting down..." << '\n';
    close(sockfd);
}
