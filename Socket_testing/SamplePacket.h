#pragma once
#include <cstdint>
#include <vector>

struct sample_packet {
    uint8_t FrameType;
    uint8_t MainUnitNum;
    uint8_t Reserved[2];
    uint32_t PacketSeqNo;
    uint16_t NumChannels;
    uint16_t NumSampleBundles;
    uint64_t FirstSampleIndex;
    uint64_t FirstSampleTime;
    // Omitting Samples due to complexity in variable-length handling
};

sample_packet deserializeSamplePacket(const std::vector<uint8_t>& buffer);
void printSamplePacket(const sample_packet& packet);
