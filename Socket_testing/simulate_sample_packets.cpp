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

std::vector<uint8_t> serializeSamplePacketData(
    uint8_t type1, uint8_t type2, uint8_t type3[2], uint32_t type4, 
    uint16_t type5, uint16_t type6, uint64_t type7, uint64_t type8, 
    std::vector<std::vector<int32_t>> &int24Data) {
    
    std::vector<uint8_t> buffer;
    // Reserve buffer space estimation (adjust based on int24Data size)
    buffer.reserve(1472);

    // Serialize fixed-size fields
    buffer.push_back(type1);
    buffer.push_back(type2);
    buffer.insert(buffer.end(), type3, type3 + 2); // Assuming type3 is an array
    // For uint32, uint16, and uint64, consider endianness
    // Here we're directly pushing back data for simplicity
    // You'll need to serialize these considering your platform's endianness
    
    // Example for uint32_t, similar approach for uint16_t and uint64_t
    for(int i = 0; i < 4; ++i) {
        buffer.push_back((type4 >> (i * 8)) & 0xFF);
    }

    for(int i = 0; i < 2; ++i) {
        buffer.push_back((type5 >> (i * 8)) & 0xFF);
    }

    for(int i = 0; i < 2; ++i) {
        buffer.push_back((type6 >> (i * 8)) & 0xFF);
    }

    for(int i = 0; i < 8; ++i) {
        buffer.push_back((type7 >> (i * 8)) & 0xFF);
    }

    for(int i = 0; i < 8; ++i) {
        buffer.push_back((type8 >> (i * 8)) & 0xFF);
    }

    // Serialize int24[N][C]
    for (auto &row : int24Data) {
        for (auto &val : row) {
            // Serialize each int32_t value as int24, adjust for your needs
            buffer.push_back((val >> 16) & 0xFF);
            buffer.push_back((val >> 8) & 0xFF);
            buffer.push_back(val & 0xFF);
        }
    }

    return buffer;
}

// Example usage
std::vector<uint8_t> generateExampleSamplePacket() {
    uint8_t FrameType = 2, MainUnitNum = 2, Reserved[2] = {3, 4};
    uint32_t PacketSeqNo = 123456;
    uint16_t NumChannels = 123, NumSampleBundles = 456;
    uint64_t FirstSampleIndex = 123456789012345, FirstSampleTime = 98765432109876;

    std::vector<std::vector<int32_t>> Samples = {
        {1234567, -1234567, 1234567, -1234567}
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

    for(int i = 0; i < 10; i++) {
        std::vector<uint8_t> data = generateExampleSamplePacket();

        sendUDP(data, "127.0.0.1", 8080);

        std::cout << "Package " << i << " sent!" << '\n';

        std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    
    return 0;
}































// std::vector<uint8_t> serializeMeasurementStartPacketData(
//     uint8_t type1, uint8_t type2, uint8_t type3[2], uint32_t type4, 
//     uint32_t type5, uint32_t type6, uint16_t type7, std::vector<std::vector<int16_t>> &type8, 
//     std::vector<std::vector<int8_t>> &type9) {
    
//     std::vector<uint8_t> buffer;
//     // Reserve buffer space estimation (adjust based on int24Data size)
//     buffer.reserve(1472);

//     // Serialize fixed-size fields
//     buffer.push_back(type1);
//     buffer.push_back(type2);
//     buffer.insert(buffer.end(), type3, type3 + 2); // Assuming type3 is an array
//     // For uint32, uint16, and uint64, consider endianness
//     // Here we're directly pushing back data for simplicity
//     // You'll need to serialize these considering your platform's endianness
    
//     // Example for uint32_t, similar approach for uint16_t and uint64_t
//     for(int i = 0; i < 4; ++i) {
//         buffer.push_back((type4 >> (i * 8)) & 0xFF);
//     }

//     // Serialize int24[N][C]
//     for (auto &row : type9) {
//         for (auto &val : row) {
//             // Serialize each int32_t value as int24, adjust for your needs
//             buffer.push_back((val >> 16) & 0xFF);
//             buffer.push_back((val >> 8) & 0xFF);
//             buffer.push_back(val & 0xFF);
//         }
//     }

//     return buffer;
// }



// Example usage
// std::vector<uint8_t> generateExampleMeasurementStartPacket() {
//     uint8_t FrameType = 1, MainUnitNum = 2, Reserved[2] = {3, 4};
//     uint32_t SamplingRateHz = 123456, SampleFormat = 234561, TriggerDefs = 345612;
//     uint16_t NumChannels = 123;
//     std::vector<std::vector<int32_t>> SourceChannels = {
//         {1234567, -1234567}, {2345678, -2345678}
//     };

//     std::vector<std::vector<int32_t>> ChannelTypes = {
//         {1234567, -1234567}, {2345678, -2345678}
//     };

// }
