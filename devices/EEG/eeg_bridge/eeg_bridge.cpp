#include "eeg_bridge.h"

// AMPLIFIER *100
// NANO TO MICRO /1000

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
    timeout.tv_sec = socket_timeout; // Timeout after 10 second
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

void EegBridge::spin(volatile std::sig_atomic_t &signal_received) {
    try {
    
    running = true;
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
        case 0x01: { // MeasurementStartPacket
            
            std::cout << "MeasurementStart package received!\n";
            measurement_start_packet packet_info;
            std::vector<uint16_t> SourceChannels;
            std::vector<uint8_t> ChannelTypes;
            
            deserializeMeasurementStartPacket_pointer(buffer, n, packet_info, SourceChannels, ChannelTypes);

            // Divide channels into data and trigger sources
            std::vector<uint16_t> data_channel_sources;
            uint16_t trigger_channel_source = -1;
            for (size_t i = 0; i < SourceChannels.size(); i++) {
                uint16_t source = SourceChannels[i];
                if(source < 60000) { 
                    data_channel_sources.push_back(source);
                } else {
                    trigger_channel_source = source; 
                }
            }

            numChannels = packet_info.NumChannels;
            
            if (trigger_channel_source != -1) numDataChannels = numChannels - 1;              // Excluding trigger channel
            else numDataChannels = numChannels;
            
            sampling_rate = packet_info.SamplingRateHz;
            lastSequenceNumber = -1;


            data_handler_samples = Eigen::MatrixXd::Zero(numDataChannels, 100);

            std::cout << "MeasurementStart package processed!\n";

            std::cout << "DataHandler reset!\n";
            eeg_bridge_status = MEASUREMENT_IN_PROGRESS;
            std::cout << "Waiting for packets..." << '\n';

            break;

        } case 0x02: { // SamplesPacket
            if (eeg_bridge_status == WAITING_MEASUREMENT_START) break;

            // Deserialize the received data into a sample_packet instance
            sample_packet packet_info;
            Eigen::VectorXi triggers_A;
            Eigen::VectorXi triggers_B;

            deserializeSamplePacketEigen_pointer(buffer, n, packet_info, data_handler_samples, triggers_A, triggers_B, (numChannels > numDataChannels));
            
            int sequenceNumber = packet_info.PacketSeqNo;

            if (numDataChannels != data_handler_samples.rows()) {
                std::cerr << "Error: numDataChannels is not equal to total rows in data_handler_samples. " << numDataChannels << ' ' << data_handler_samples.rows() << '\n';
                break;
            }

            Eigen::MatrixXd data_samples = ((data_handler_samples * DC_MODE_SCALE) / NANO_TO_MICRO_CONVERSION);

            // Ensure column access is valid
            if (packet_info.NumSampleBundles > data_samples.cols()) {
                std::cerr << "Error: Requested more sample bundles than available columns in data_samples." << '\n';
                break;
            }

            for (int i = 0; i < packet_info.NumSampleBundles; i++) {
                if (i >= data_samples.cols()) {
                    std::cerr << "Error: Column index " << i << " is out of range for data_samples with columns " << data_samples.cols() << '\n';
                    break;
                }
            }

            // Check for dropped packets
            if (lastSequenceNumber != -1 && sequenceNumber != (lastSequenceNumber + 1)) {
                std::cerr << "Packet loss detected. Expected sequence: " << (lastSequenceNumber + 1) << ", but received: " << sequenceNumber << '\n';
                // TODO: Handle packet loss
            }

            lastSequenceNumber = sequenceNumber; // Update the latest sequence number

            // Debug output to confirm data integrity
            // std::cout << "Package " << sequenceNumber << " received!\n";
            // std::cout << "Channels: " << data_samples.rows() << ", " << data_samples.cols() << '\n';
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

    running = false;
    std::cout << "Shutting down..." << '\n';
    close(sockfd);

    } catch (const std::exception& e) {
        std::cerr << "Eeg_bridge exception: " << e.what() << '\n';
        std::cerr << boost::stacktrace::stacktrace();
    }
}
