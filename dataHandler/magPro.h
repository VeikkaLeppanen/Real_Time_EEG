#ifndef MAGPRO_H
#define MAGPRO_H

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <fstream>
#include <boost/asio.hpp>

class magPro {
public:
    magPro()
        : serial(io)
    {
    }

    int connectTriggerPort();
    void trig();
    std::vector<uint8_t> create_trig_cmd_byte_str();
    std::vector<uint8_t> create_enable_cmd_byte_str(bool enable);
    std::vector<uint8_t> create_amplitude_cmd_byte_str(uint8_t amplitudeA_value, uint8_t amplitudeB_value);
    uint8_t crc8(const std::vector<uint8_t>& data);
    void set_enable(bool status);
    void set_amplitude(int amplitude);

private:
    boost::asio::io_service io;
    boost::asio::serial_port serial;
};

#endif // MAGPRO_H