#include "samplePacket.h"

void deserializeSamplePacketEigen_pointer(const uint8_t *buffer, size_t size, sample_packet &packet, Eigen::MatrixXd &packet_handler_buffer) {
    size_t offset = 0;

    if (buffer == nullptr || size < sizeof(sample_packet)) {
        throw std::runtime_error("Invalid buffer or size.");
    }

    packet.FrameType = buffer[offset++];
    packet.MainUnitNum = buffer[offset++];
    packet.Reserved[0] = buffer[offset++];
    packet.Reserved[1] = buffer[offset++];

    memcpy(&packet.PacketSeqNo, buffer + offset, sizeof(uint32_t));
    packet.PacketSeqNo = ntohl(packet.PacketSeqNo);
    offset += sizeof(uint32_t);

    memcpy(&packet.NumChannels, buffer + offset, sizeof(uint16_t));
    packet.NumChannels = ntohs(packet.NumChannels);
    offset += sizeof(uint16_t);

    memcpy(&packet.NumSampleBundles, buffer + offset, sizeof(uint16_t));
    packet.NumSampleBundles = ntohs(packet.NumSampleBundles);
    offset += sizeof(uint16_t);

    memcpy(&packet.FirstSampleIndex, buffer + offset, sizeof(uint64_t));
    packet.FirstSampleIndex = ntohll(packet.FirstSampleIndex);
    offset += sizeof(uint64_t);

    memcpy(&packet.FirstSampleTime, buffer + offset, sizeof(uint64_t));
    packet.FirstSampleTime = ntohll(packet.FirstSampleTime);
    offset += sizeof(uint64_t);

    // Initialize an Eigen matrix to store samples. Rows correspond to channels, columns to sample bundles.

    for (uint16_t bundle = 0; bundle < packet.NumSampleBundles; ++bundle) {
        for (uint16_t channel = 0; channel < packet.NumChannels; ++channel) {
            if (offset + 3 > size) {
                throw std::runtime_error("Buffer overrun while reading samples.");
            }
            // Extracting and converting 24-bit signed integer to double
            int32_t sample = (static_cast<int32_t>(buffer[offset]) << 24) |
                             (static_cast<int32_t>(buffer[offset + 1]) << 16) |
                             (static_cast<int32_t>(buffer[offset + 2]) << 8);
            sample >>= 8; // Sign extension for 24-bit to 32-bit
            packet_handler_buffer(channel, bundle) = static_cast<double>(sample);
            offset += 3;
        }
    }

    return;
}


std::vector<std::vector<double>> deserializeSamplePacket_pointer(const uint8_t *buffer, size_t size, sample_packet &packet) {
    std::vector<std::vector<double>> data;
    size_t offset = 0;

    if (buffer == nullptr || size < sizeof(sample_packet)) {
        // Handle error: invalid buffer or size
        throw std::runtime_error("Invalid buffer or size.");
    }

    // Direct copy for single byte fields
    packet.FrameType = buffer[offset++];
    packet.MainUnitNum = buffer[offset++];
    packet.Reserved[0] = buffer[offset++];
    packet.Reserved[1] = buffer[offset++];

    // Assuming big-endian network byte order for multi-byte fields
    memcpy(&packet.PacketSeqNo, buffer + offset, sizeof(uint32_t));
    packet.PacketSeqNo = ntohl(packet.PacketSeqNo);
    offset += sizeof(uint32_t);

    memcpy(&packet.NumChannels, buffer + offset, sizeof(uint16_t));
    packet.NumChannels = ntohs(packet.NumChannels);
    offset += sizeof(uint16_t);

    memcpy(&packet.NumSampleBundles, buffer + offset, sizeof(uint16_t));
    packet.NumSampleBundles = ntohs(packet.NumSampleBundles);
    offset += sizeof(uint16_t);

    memcpy(&packet.FirstSampleIndex, buffer + offset, sizeof(uint64_t));
    packet.FirstSampleIndex = ntohll(packet.FirstSampleIndex);
    offset += sizeof(uint64_t);

    memcpy(&packet.FirstSampleTime, buffer + offset, sizeof(uint64_t));
    packet.FirstSampleTime = ntohll(packet.FirstSampleTime);
    offset += sizeof(uint64_t);

    // Initialize data structure to store samples
    data.resize(packet.NumChannels);

    // Handle Samples based on NumChannels and NumSampleBundles, if applicable
    for (uint16_t bundle = 0; bundle < packet.NumSampleBundles; ++bundle) {
        for (uint16_t channel = 0; channel < packet.NumChannels; ++channel) {
            if (offset + 3 > size) { // Ensure there are 3 bytes available for each sample
                throw std::runtime_error("Buffer overrun while reading samples.");
            }
            // Extracting and converting 24-bit signed integer to double
            int32_t sample = (buffer[offset] << 24) | (buffer[offset + 1] << 16) | (buffer[offset + 2] << 8);
            sample >>= 8; // Sign extension for 24-bit to 32-bit
            double sampleDouble = static_cast<double>(sample);
            data[channel].push_back(sampleDouble);
            offset += 3; // Move past the 24-bit sample
        }
    }

    return data;
}

void printSamplePacket(const sample_packet& packet) {
    std::cout << "FrameType: " << static_cast<int>(packet.FrameType) << std::endl;
    std::cout << "MainUnitNum: " << static_cast<int>(packet.MainUnitNum) << std::endl;
    std::cout << "Reserved: " << static_cast<int>(packet.Reserved[0]) 
              << ", " << static_cast<int>(packet.Reserved[1]) << std::endl;
    std::cout << "PacketSeqNo: " << packet.PacketSeqNo << std::endl;
    std::cout << "NumChannels: " << packet.NumChannels << std::endl;
    std::cout << "NumSampleBundles: " << packet.NumSampleBundles << std::endl;
    std::cout << "FirstSampleIndex: " << packet.FirstSampleIndex << std::endl;
    std::cout << "FirstSampleTime: " << packet.FirstSampleTime << std::endl;
    std::cout << std::endl;
}
