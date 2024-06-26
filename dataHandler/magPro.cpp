#include "magPro.h"

int magPro::connectTriggerPort() {

    // Check if the port exists
    std::string port = "/dev/ttyUSB0";
    if (!std::ifstream(port)) {
        std::cerr << "Error: Port " << port << " not found." << std::endl;
        return 1; // Return error if port does not exist
    }

    // Create a serial port object
    serial.open(port);  // Specify the serial port to connect to

    try {
        // Set the baud rate
        serial.set_option(boost::asio::serial_port_base::baud_rate(38400));

        // Set parity (none)
        serial.set_option(boost::asio::serial_port_base::parity(boost::asio::serial_port_base::parity::none));

        // Set stop bits (one)
        serial.set_option(boost::asio::serial_port_base::stop_bits(boost::asio::serial_port_base::stop_bits::one));

        // Set character size (8 bits)
        serial.set_option(boost::asio::serial_port_base::character_size(8));

        std::cout << "Serial port configured and opened." << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));

    } catch (const boost::system::system_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;  // Return success
}

void magPro::trig() {

    auto cmd = create_trig_cmd_byte_str();
    boost::asio::write(serial, boost::asio::buffer(cmd));
    
}

void magPro::set_enable(bool status) {
    auto cmd_enable = create_enable_cmd_byte_str(status);
    boost::asio::write(serial, boost::asio::buffer(cmd_enable));
    std::this_thread::sleep_for(std::chrono::seconds(1));
}

void magPro::set_amplitude(int amplitude) {
    if (amplitude <= 0 || amplitude >= 100) {
        std::cerr << "Error amplitude must be between 0 and 100" << std::endl;
        return;
    } else {
        auto cmd_amplitude = create_amplitude_cmd_byte_str(amplitude, 0);
        boost::asio::write(serial, boost::asio::buffer(cmd_amplitude));
    }
    std::this_thread::sleep_for(std::chrono::seconds(1));
}

std::vector<uint8_t> magPro::create_trig_cmd_byte_str() {
    std::vector<uint8_t> cmd = {0x03, 0x01, 0x00};
    uint8_t crc = crc8(cmd);
    std::vector<uint8_t> fullCmd = {0xfe, static_cast<uint8_t>(cmd.size())};
    fullCmd.reserve(4 + cmd.size()); // 4 extra for {0xfe, size, crc, 0xff}
    fullCmd.insert(fullCmd.end(), cmd.begin(), cmd.end());
    fullCmd.push_back(crc);
    fullCmd.push_back(0xff);
    return fullCmd;
}

std::vector<uint8_t> magPro::create_enable_cmd_byte_str(bool enable) {
    std::vector<uint8_t> cmd;

    if (enable) {
        cmd = {0x02, 0x01, 0x00};
    } else {
        cmd = {0x02, 0x00, 0x00};
    }

    uint8_t crc = crc8(cmd);
    std::vector<uint8_t> result = {0xfe, static_cast<uint8_t>(cmd.size())};
    result.reserve(4 + cmd.size());
    result.insert(result.end(), cmd.begin(), cmd.end());
    result.push_back(crc);
    result.push_back(0xff);

    return result;
}

std::vector<uint8_t> magPro::create_amplitude_cmd_byte_str(uint8_t amplitudeA_value, uint8_t amplitudeB_value) {
    std::vector<uint8_t> cmd = {0x01, amplitudeA_value, amplitudeB_value};

    uint8_t crc = crc8(cmd);
    std::vector<uint8_t> result = {0xfe, static_cast<uint8_t>(cmd.size())};
    result.reserve(4 + cmd.size());
    result.insert(result.end(), cmd.begin(), cmd.end());
    result.push_back(crc);
    result.push_back(0xff);

    return result;
}

uint8_t magPro::crc8(const std::vector<uint8_t>& bytes) {
    uint8_t crc = 0;
    for (auto next : bytes) {
        for (int j = 0; j < 8; j++) {
            if ((next ^ crc) & 0x01) {
                crc ^= 0x18;
                crc = (crc & 0xff) >> 1;
                crc |= 0x80;
            } else {
                crc >>= 1;
            }
            next >>= 1;
        }
    }
    return crc;
}