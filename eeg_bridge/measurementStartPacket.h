#pragma once
#include "networkUtils.h"
#include <iostream>
#include <arpa/inet.h>
#include <cstdint>
#include <cstring>

struct __attribute__((packed)) measurement_start_packet {
    uint8_t FrameType;
    uint8_t MainUnitNum;
    uint8_t Reserved[2];
    uint32_t SamplingRateHz;
    uint32_t SampleFormat;
    uint32_t TriggerDefs;
    uint16_t NumChannels;
};

void deserializeMeasurementStartPacket_pointer(const uint8_t *buffer, size_t size, measurement_start_packet &packet);
void printMeasurementStartPacket(const measurement_start_packet &packet);
