#include "measurementStartPacket.h"

void deserializeMeasurementStartPacket_pointer(const uint8_t *buffer, size_t size, measurement_start_packet &packet, std::vector<uint16_t> &SourceChannels, std::vector<uint8_t> &ChannelTypes) {
    size_t offset = 0;

    if (buffer == nullptr || size < sizeof(measurement_start_packet)) {
        // Handle error: invalid buffer or size
        throw std::runtime_error("Invalid buffer or size.");
    }

    // Direct copy for single byte fields
    packet.FrameType = buffer[offset++];
    packet.MainUnitNum = buffer[offset++];
    packet.Reserved[0] = buffer[offset++];
    packet.Reserved[1] = buffer[offset++];

    // Assuming big-endian network byte order for multi-byte fields
    memcpy(&packet.SamplingRateHz, buffer + offset, sizeof(uint32_t));
    packet.SamplingRateHz = ntohl(packet.SamplingRateHz);
    offset += sizeof(uint32_t);

    memcpy(&packet.SampleFormat, buffer + offset, sizeof(uint32_t));
    packet.SampleFormat = ntohl(packet.SampleFormat);
    offset += sizeof(uint32_t);

    memcpy(&packet.TriggerDefs, buffer + offset, sizeof(uint32_t));
    packet.TriggerDefs = ntohl(packet.TriggerDefs);
    offset += sizeof(uint32_t);

    memcpy(&packet.NumChannels, buffer + offset, sizeof(uint16_t));
    packet.NumChannels = ntohs(packet.NumChannels);
    offset += sizeof(uint16_t);
    
    // Source channels and channel types can be commented in if needed
    
    // Assuming SourceChannels are next in the buffer
    SourceChannels.resize(packet.NumChannels);
    for (uint16_t channel = 0; channel < packet.NumChannels; ++channel) {
        if (offset + sizeof(uint16_t) > size) throw std::runtime_error("Buffer overrun while reading SourceChannels.");
        uint16_t sourceChannel;
        memcpy(&sourceChannel, buffer + offset, sizeof(uint16_t));
        SourceChannels[channel] = ntohs(sourceChannel);
        offset += sizeof(uint16_t);
    }
    
    // Assuming ChannelTypes follow SourceChannels in the buffer
    ChannelTypes.resize(packet.NumChannels);
    for (uint16_t channel = 0; channel < packet.NumChannels; ++channel) {
        if (offset + sizeof(uint8_t) > size) throw std::runtime_error("Buffer overrun while reading ChannelTypes.");
        ChannelTypes[channel] = buffer[offset++];
    }
}

void printMeasurementStartPacket(const measurement_start_packet &packet) {
    std::cout << "FrameType: " << static_cast<int>(packet.FrameType) << std::endl;
    std::cout << "MainUnitNum: " << static_cast<int>(packet.MainUnitNum) << std::endl;
    std::cout << "Reserved: " << static_cast<int>(packet.Reserved[0]) 
              << ", " << static_cast<int>(packet.Reserved[1]) << std::endl;
    std::cout << "SamplingRateHz: " << packet.SamplingRateHz << std::endl;
    std::cout << "SampleFormat: " << packet.SampleFormat << std::endl;
    std::cout << "PacketSeqNo: " << packet.TriggerDefs << std::endl;
    std::cout << "NumChannels: " << packet.NumChannels << std::endl;
    std::cout << std::endl;
}
