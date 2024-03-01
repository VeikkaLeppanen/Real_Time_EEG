#include <iostream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <chrono>
#include <thread>

/*
This can be used to send packets to the port used by the main program to simulate the retrieval of bittium packets.
Current available packages:
    MeasurementStartPackage
    SamplesPackage

TODO:
    MeasurementEndPackage
*/

// Target port
#define PORT 8080


std::vector<uint8_t> serializeMeasurementStartPacketData(
    uint8_t FrameType,
    uint8_t MainUnitNum,
    uint8_t Reserved[2],
    uint32_t SamplingRateHz,
    uint32_t SampleFormat,
    uint32_t TriggerDefs,
    uint16_t NumChannels) {
    
    std::vector<uint8_t> buffer;
    
    // Assuming the buffer should have a basic capacity. Adjust as needed.
    buffer.reserve(1472); // Reserve an arbitrary size for simplicity

    // Serialize fixed-size fields in big-endian order
    buffer.push_back(FrameType);
    buffer.push_back(MainUnitNum);
    buffer.insert(buffer.end(), Reserved, Reserved + 2); // Assuming Reserved is an array
    
    // Serialize multi-byte integers in big-endian order
    buffer.push_back((SamplingRateHz >> 24) & 0xFF);
    buffer.push_back((SamplingRateHz >> 16) & 0xFF);
    buffer.push_back((SamplingRateHz >> 8) & 0xFF);
    buffer.push_back(SamplingRateHz & 0xFF);

    buffer.push_back((SampleFormat >> 24) & 0xFF);
    buffer.push_back((SampleFormat >> 16) & 0xFF);
    buffer.push_back((SampleFormat >> 8) & 0xFF);
    buffer.push_back(SampleFormat & 0xFF);

    buffer.push_back((TriggerDefs >> 24) & 0xFF);
    buffer.push_back((TriggerDefs >> 16) & 0xFF);
    buffer.push_back((TriggerDefs >> 8) & 0xFF);
    buffer.push_back(TriggerDefs & 0xFF);

    buffer.push_back((NumChannels >> 8) & 0xFF);
    buffer.push_back(NumChannels & 0xFF);

    return buffer;
}

std::vector<uint8_t> generateExampleMeasurementStartPacket() {
    uint8_t FrameType = 1; // Assuming 1 indicates a Measurement Start packet
    uint8_t MainUnitNum = 2;
    uint8_t Reserved[2] = {0, 0}; // Fill with appropriate values if needed
    uint32_t SamplingRateHz = 5000; // Example: 5000 Hz
    uint32_t SampleFormat = 1; // Example format
    uint32_t TriggerDefs = 0; // Example trigger definitions
    uint16_t NumChannels = 4; // Number of channels

    return serializeMeasurementStartPacketData(
        FrameType, MainUnitNum, Reserved, SamplingRateHz,
        SampleFormat, TriggerDefs, NumChannels);
}

std::vector<uint8_t> serializeSamplePacketData(
    uint8_t FrameType, uint8_t MainUnitNum, uint8_t Reserved[2], uint32_t PacketSeqNo, 
    uint16_t NumChannels, uint16_t NumSampleBundles, uint64_t FirstSampleIndex, uint64_t FirstSampleTime, 
    std::vector<std::vector<int32_t>> &int24Data) {
    
    std::vector<uint8_t> buffer;
    buffer.reserve(1472); // Adjust based on int24Data size

    // Serialize fixed-size fields in big-endian order
    buffer.push_back(FrameType);
    buffer.push_back(MainUnitNum);
    buffer.insert(buffer.end(), Reserved, Reserved + 2);
    
    // Serialize multi-byte integers in big-endian order
    buffer.push_back((PacketSeqNo >> 24) & 0xFF);
    buffer.push_back((PacketSeqNo >> 16) & 0xFF);
    buffer.push_back((PacketSeqNo >> 8) & 0xFF);
    buffer.push_back(PacketSeqNo & 0xFF);

    buffer.push_back((NumChannels >> 8) & 0xFF);
    buffer.push_back(NumChannels & 0xFF);

    buffer.push_back((NumSampleBundles >> 8) & 0xFF);
    buffer.push_back(NumSampleBundles & 0xFF);

    for(int i = 7; i >= 0; --i) {
        buffer.push_back((FirstSampleIndex >> (i * 8)) & 0xFF);
    }

    for(int i = 7; i >= 0; --i) {
        buffer.push_back((FirstSampleTime >> (i * 8)) & 0xFF);
    }

    // Serialize int24[N][C] in big-endian order
    for (auto &row : int24Data) {
        for (auto &val : row) {
            // Assuming val represents a signed 24-bit integer correctly, including handling negative values
            buffer.push_back((val >> 16) & 0xFF); // MSB
            buffer.push_back((val >> 8) & 0xFF);
            buffer.push_back(val & 0xFF); // LSB
        }
    }

    return buffer;
}

// Example usage
std::vector<uint8_t> generateExampleSamplePacket() {
    uint8_t FrameType = 2, MainUnitNum = 2, Reserved[2] = {3, 4};
    uint32_t PacketSeqNo = 123456;
    uint16_t NumChannels = 4, NumSampleBundles = 2;
    uint64_t FirstSampleIndex = 123456789012345, FirstSampleTime = 98765432109876;

    std::vector<std::vector<int32_t>> Samples = {
        {1234567, -1234567}, {1234567, -1234567}, {1234567, -1234567}, {1234567, -1234567}
    };

    return serializeSamplePacketData(FrameType, MainUnitNum, Reserved, PacketSeqNo, NumChannels, NumSampleBundles, FirstSampleIndex, FirstSampleTime, Samples);
}

void sendUDP(const std::vector<uint8_t> &data, const std::string &address, int port) {
    int sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0) {
        std::cerr << "Error opening socket" << std::endl;
        return;
    }

    sockaddr_in servaddr;
    memset(&servaddr, 0, sizeof(servaddr));
    servaddr.sin_family = AF_INET;
    servaddr.sin_port = htons(port);
    servaddr.sin_addr.s_addr = inet_addr(address.c_str());

    sendto(sockfd, data.data(), data.size(), 0, (const sockaddr *)&servaddr, sizeof(servaddr));
    close(sockfd);
}

int main() {

    std::vector<uint8_t> MSdata = generateExampleMeasurementStartPacket();

    sendUDP(MSdata, "127.0.0.1", PORT);

    std::cout << "MeasurementStartPackage sent!" << '\n';

    std::this_thread::sleep_for(std::chrono::seconds(1));

    for(int i = 0; i < 3; i++) {
        std::vector<uint8_t> sample_data = generateExampleSamplePacket();

        sendUDP(sample_data, "127.0.0.1", PORT);

        std::cout << "Package " << i << " sent!" << '\n';

        std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    
    return 0;
}
