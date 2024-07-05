#ifndef MAGPRO_H
#define MAGPRO_H

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <chrono>
#include <fstream>
#include <boost/asio.hpp>

class magPro {
public:
    magPro()
        : serial(io)
    { 
        latest_trigger_time = std::chrono::system_clock::now();
    }

    int connectTriggerPort();
    void trig();
    std::vector<uint8_t> create_trig_cmd_byte_str();
    std::vector<uint8_t> create_enable_cmd_byte_str(bool enable);
    std::vector<uint8_t> create_amplitude_cmd_byte_str(uint8_t amplitudeA_value, uint8_t amplitudeB_value);
    uint8_t crc8(const std::vector<uint8_t>& data);
    void set_enable(bool status);
    void set_amplitude(int amplitude);
    void setTriggerTimeLimit(int value) { time_limit = std::max(min_time_limit, std::min(max_time_limit, value)); }
    bool getTriggerTimeLimit() { return time_limit; }

private:
    boost::asio::io_service io;
    boost::asio::serial_port serial;

    // Time limit in milliseconds
    int time_limit;
    const int min_time_limit = 100;
    const int max_time_limit = 10000;
    std::chrono::time_point<std::chrono::system_clock> latest_trigger_time;
};

#endif // MAGPRO_H