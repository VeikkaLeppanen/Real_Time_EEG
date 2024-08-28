#pragma once
#include "networkUtils.h"
#include <iostream>
#include <arpa/inet.h>
#include <stdexcept>
#include <cstdint>
#include <cstring>
#include <Eigen/Dense>

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

void deserializeSamplePacketEigen_pointer(const uint8_t *buffer, size_t size, sample_packet &packet, Eigen::MatrixXd &packet_handler_buffer, Eigen::VectorXi &triggers_A, Eigen::VectorXi &triggers_B, bool contains_trigger_channel);
std::vector<std::vector<double>> deserializeSamplePacket_pointer(const uint8_t *buffer, size_t size, sample_packet &packet);
void printSamplePacket(const sample_packet& packet);
