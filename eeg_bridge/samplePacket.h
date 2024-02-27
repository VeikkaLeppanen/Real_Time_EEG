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
};

sample_packet deserializeSamplePacket_vector(const std::vector<uint8_t> &buffer);
std::vector<std::vector<double>> deserializeSamplePacket_pointer(const uint8_t *buffer, size_t size, sample_packet &packet);
void printSamplePacket(const sample_packet& packet);
