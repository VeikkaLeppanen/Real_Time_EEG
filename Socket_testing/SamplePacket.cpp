#include "SamplePacket.h"
#include "NetworkUtils.h"
#include <iostream>
#include <arpa/inet.h>

sample_packet deserializeSamplePacket(const std::vector<uint8_t>& buffer) {
    sample_packet packet;
    size_t offset = 0;

    // Direct copy for single byte fields
    packet.FrameType = buffer[offset++];
    packet.MainUnitNum = buffer[offset++];
    packet.Reserved[0] = buffer[offset++];
    packet.Reserved[1] = buffer[offset++];

    // Assuming big-endian network byte order for multi-byte fields
    packet.PacketSeqNo = ntohl(*reinterpret_cast<const uint32_t*>(&buffer[offset]));
    offset += sizeof(uint32_t);

    packet.NumChannels = ntohs(*reinterpret_cast<const uint16_t*>(&buffer[offset]));
    offset += sizeof(uint16_t);

    packet.NumSampleBundles = ntohs(*reinterpret_cast<const uint16_t*>(&buffer[offset]));
    offset += sizeof(uint16_t);

    packet.FirstSampleIndex = ntohll(*reinterpret_cast<const uint64_t*>(&buffer[offset]));
    offset += sizeof(uint64_t);

    packet.FirstSampleTime = ntohll(*reinterpret_cast<const uint64_t*>(&buffer[offset]));
    offset += sizeof(uint64_t);

    // Handle Samples based on NumChannels and NumSampleBundles, if applicable

    return packet;
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
