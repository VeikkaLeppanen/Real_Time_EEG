#include "eeg_bridge.h"

const int PORT = 50000;

// AMPLIFIER *100
// NANO TO MICRO /1000
const uint8_t DC_MODE_SCALE = 100;
const uint16_t NANO_TO_MICRO_CONVERSION = 1000;

void EegBridge::bind_socket() {

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

    struct timeval timeout;
    timeout.tv_sec = 10; // Timeout after 10 second
    timeout.tv_usec = 0;

    if (setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, (char *)&timeout, sizeof(timeout)) < 0) {
        std::cerr << "Failed to set socket timeout." << std::endl;
    }

  
    // Set socket buffer size to 10 MB
    int buffer_size = 1024 * 1024 * 10;
    if (setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, &buffer_size, sizeof(buffer_size)) == -1) {
        std::cerr << "Failed to set socket receive buffer size." << std::endl;
    }

    len = sizeof(cliaddr);
}

void EegBridge::spin(dataHandler &handler, volatile std::sig_atomic_t &signal_received) {
    std::cout << "Waiting for measurement start..." << '\n';
    eeg_bridge_status = WAITING_MEASUREMENT_START;
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

            if (eeg_bridge_status == WAITING_MEASUREMENT_START) break;
        
            // Deserialize the received data into a sample_packet instance
            sample_packet packet_info;
            deserializeSamplePacketEigen_pointer(buffer, n, packet_info, data_handler_samples);

            Eigen::VectorXd triggers = data_handler_samples.row(data_handler_samples.rows() - 1);

            Eigen::MatrixXd data_samples = (data_handler_samples.topRows(data_handler_samples.rows() - 1) * DC_MODE_SCALE) / NANO_TO_MICRO_CONVERSION;

            for (int i = 0; i < packet_info.NumSampleBundles; i++) {
                handler.addData(data_samples.col(i), static_cast<double>(packet_info.FirstSampleTime), static_cast<int>(triggers(i)));
            }

            int sequenceNumber = packet_info.PacketSeqNo;

            // Check for dropped packets
            if (lastSequenceNumber != -1 && sequenceNumber != (lastSequenceNumber + 1)) {
                std::cerr << "Packet loss detected. Expected sequence: " << (lastSequenceNumber + 1) << ", but received: " << sequenceNumber << '\n';
                // Handle packet loss (e.g., by requesting retransmission if supported)
            }

            lastSequenceNumber = sequenceNumber; // Update the last received sequence number

            std::cout << "Package " << sequenceNumber << " received!" << '\n';

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
            numDataChannels = numChannels - 1;              // Excluding trigger channel
            sampling_rate = packet_info.SamplingRateHz;
            lastSequenceNumber = -1;

            // TODO: Initialize data_handler_samples
            data_handler_samples = Eigen::MatrixXd::Zero(numChannels, 10);


            std::cout << "MeasurementStart package processed!\n";

            handler.reset_handler(numDataChannels, sampling_rate);
            std::cout << "DataHandler reset!\n";
            eeg_bridge_status = MEASUREMENT_IN_PROGRESS;
            std::cout << "Waiting for packets..." << '\n';

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
