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

    // Check if the trigger has been triggered recently
    if (std::chrono::system_clock::now() - latest_trigger_time < std::chrono::duration<double>(time_limit)) { return; }
    
    auto cmd = create_trig_cmd_byte_str();
    boost::asio::write(serial, boost::asio::buffer(cmd));

    sleep(sleep_time_trig);
}

void magPro::set_enable(bool status) {
    auto cmd_enable = create_enable_cmd_byte_str(status);
    boost::asio::write(serial, boost::asio::buffer(cmd_enable));

    sleep(sleep_time_set_enable);
}

void magPro::set_amplitude(int amplitude) {
    if (amplitude <= 0 || amplitude >= 100) {
        std::cerr << "Error amplitude must be between 0 and 100" << std::endl;
        return;
    } else {
        auto cmd_amplitude = create_amplitude_cmd_byte_str(amplitude, 0);
        boost::asio::write(serial, boost::asio::buffer(cmd_amplitude));
    }

    sleep(sleep_time_set_amp);
}






// Data receiving functions
int magPro::wait_for_next_package_from_G3() {
        std::string received_cmd;
        bool received_end_of_cmd = false;
        int received_cmd_length;
        while (!received_end_of_cmd) {
            unsigned char read_byte;
            read_from_serial(&read_byte, 1);
            print_debug("Received: " + read_byte);
            if (read_byte == 0xFE) { // start of command registered
                print_debug("Received package start 0xFE");
                unsigned char received_cmd_length_char;
                read_from_serial(&received_cmd_length_char, 1);
                received_cmd_length = static_cast<int>(received_cmd_length_char);
                print_debug("Received cmd_length=" + std::to_string(received_cmd_length));
                received_cmd.resize(received_cmd_length);
                read_from_serial(reinterpret_cast<unsigned char*>(&received_cmd[0]), received_cmd_length);
                print_debug("Received command: " + received_cmd + ": " + cmd_type_name(static_cast<CmdType>(received_cmd[0])));

                unsigned char received_crc_char;
                read_from_serial(&received_crc_char, 1);
                int received_crc = static_cast<int>(received_crc_char);
                print_debug("Received crc: " + std::to_string(received_crc) + " is of type: " + typeid(received_crc).name());

                std::vector<uint8_t> cmd_vector(received_cmd.begin(), received_cmd.end());
                int calculated_crc = crc8(cmd_vector);
                print_debug("Calculated crc: " + std::to_string(calculated_crc) + " is of type: " + typeid(calculated_crc).name());

                if (received_crc != calculated_crc) {
                    std::cerr << "crc error detected !" << std::endl;
                    return -1;
                }

                unsigned char read_end_flag;
                read_from_serial(&read_end_flag, 1);
                if (read_end_flag == 0xFF) {
                    print_debug("Received package end 0xFF");
                    received_end_of_cmd = true;
                } else {
                    std::cerr << "Error! Unknown package received. Got " << read_end_flag << " as read_end_flag" << std::endl;
                    return -1;
                }
            }
        }
        print_debug("Received command: " + received_cmd + " end");

        // General debug of CMD byte for byte
        for (size_t i = 0; i < received_cmd.size(); ++i) {
            print_debug("received_cmd[" + std::to_string(i) + "]:" + std::to_string(static_cast<int>(received_cmd[i])));
        }

        // Handling each package individually
        switch (received_cmd_length) {
            case 4:
                return handle_cmd_length_4(received_cmd);
            case 10:
                return handle_cmd_length_10(received_cmd);
            case 13:
                return handle_cmd_length_mep(received_cmd, received_cmd_length);
            case 17:
                return handle_cmd_length_mep(received_cmd, received_cmd_length);
            default:
                std::cerr << "Unsupported received_cmd_length=" + std::to_string(received_cmd_length) << std::endl;
                return -1;
        }
    }


int magPro::handle_cmd_length_4(const std::string& received_cmd) {
    switch (static_cast<CmdType>(received_cmd[0])) {
        case CmdType::AMPLITUDE:
            G3_A_amp_percent = received_cmd[1];
            G3_B_amp_percent = received_cmd[2];
            G3_amp_timestamp_ms = get_epochtime_ms();
            print_debug(std::to_string(get_epochtime_ms()) + ": Value A Amplitude in % = " + std::to_string(received_cmd[1]));
            print_debug("Value B Amplitude in % = " + std::to_string(received_cmd[2]));
            cmd_length_4__evaluate_byte_6(received_cmd[3]);
            return static_cast<int>(CmdType::AMPLITUDE);
        case CmdType::DI_DT:
            print_debug(get_timestamp_str() + ": A di/dt=" + std::to_string(received_cmd[1]) + ", B di/dt=" + std::to_string(received_cmd[2]));
            G3_A_di_dt = received_cmd[1];
            G3_B_di_dt = received_cmd[2];
            G3_di_dt_timestamp_ms = get_epochtime_ms();
            cmd_length_4__evaluate_byte_6(received_cmd[3]);
            print_formatted_data_if_all_available();
            return static_cast<int>(CmdType::DI_DT);
        case CmdType::TEMPERATURE:
            std::cout << "Ignoring received temperature package" << std::endl;
            return static_cast<int>(CmdType::TEMPERATURE);
        case CmdType::PAGE_TRAIN_RUNNING_STATUS:
            G3_current_page = received_cmd[1];
            print_debug("Stored G3_current_page=" + std::to_string(G3_current_page));
            G3_train_running = received_cmd[2];
            print_debug("Stored G3_train_running=" + std::to_string(G3_train_running));
            return static_cast<int>(CmdType::PAGE_TRAIN_RUNNING_STATUS);
        default:
            std::cerr << "Unsupported type received_cmd[0]=" + std::to_string(received_cmd[0]) << std::endl;
            return -1; // Use -1 to indicate error
    }
}

int magPro::handle_cmd_length_10(const std::string& received_cmd) {
    if (static_cast<CmdType>(received_cmd[0]) == CmdType::MODE) {
        print_debug("Matched received_cmd[0] == CmdType::MODE");
        cmd_length_10__mode_9__evaluate_byte_5(received_cmd[2]);
        store_mode(received_cmd[3]);
        store_current_direction(received_cmd[4]);
        store_waveform(received_cmd[5]);
        int burst_pulses_index = received_cmd[6];
        if (burst_pulses_index >= 0 && burst_pulses_index <= 3) {
            store_current_burst_pulses(5 - burst_pulses_index);
        } else {
            throw std::out_of_range("burst_pulses_index=" + std::to_string(burst_pulses_index) + " out of range 0-3");
        }
        store_current_ipi_index(received_cmd.substr(7, 2));
        store_current_ba_ratio(received_cmd[9]);
        print_debug("About to return CmdType::MODE");
        return static_cast<int>(CmdType::MODE);
    }
    return -1; // Use -1 to indicate error
}

int magPro::handle_cmd_length_mep(const std::string& received_cmd, int length) {
    if (static_cast<CmdType>(received_cmd[0]) == CmdType::MEP_MIN_MAX) {
        print_debug("Matched received_cmd[0] == CmdType::MEP_MIN_MAX for ExtendedMEPData=" + std::to_string(length == 13 ? 0 : 1));
        store_mep_values(received_cmd);
        return static_cast<int>(CmdType::MEP_MIN_MAX);
    }
    return -1; // Use -1 to indicate error
}

long long magPro::get_epochtime_ms() {
    auto now = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
    return duration.count();
}

std::string magPro::get_timestamp_str() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    return oss.str();
}

void magPro::cmd_length_4__evaluate_byte_6(char byte) {
    int bit_0_and_1 = byte & 3;            // mode
    int bit_2_and_3 = (byte >> 2) & 3;     // waveform
    int bit_4 = (byte >> 4) & 1;           // enabled
    int bit_5_to_7 = (byte >> 5) & 7;      // model
    
    store_mode(bit_0_and_1);
    store_waveform(bit_2_and_3);
    store_current_enabled(bit_4);
    store_model(bit_5_to_7);
}

void magPro::cmd_length_10__mode_9__evaluate_byte_5(char byte) {
    int bit_0_to_2 = byte & 7;             // model
    int bit_3_to_5 = (byte >> 3) & 7;      // Max possible pps for given model
    int bit_6 = (byte >> 6) & 1;           // enabled
    
    store_model(bit_0_to_2);
    
    if (bit_0_to_2 == 0 || bit_0_to_2 == 4) {  // R30
        switch (bit_3_to_5) {
            case 0:
                store_max_pps(30);
                break;
            case 1:
                store_max_pps(20);
                break;
            case 2:
                store_max_pps(60);
                break;
            case 3:
                store_max_pps(80);
                break;
            default:
                std::cerr << "Unsupported ppm reference in bit_3_to_5=" << bit_3_to_5 << std::endl;
        }
    } else if (bit_0_to_2 == 1 || bit_0_to_2 == 3) {  // X100
        if (bit_3_to_5 == 3) {
            store_max_pps(250);
        } else {
            store_max_pps(100);  // Default if not specified
        }
    } else {
        std::cerr << "Unsupported model bit_0_to_2=" << bit_0_to_2 << std::endl;
    }

    store_supports_biphasic_burst(bit_6);
}

void magPro::store_mode(char byte) {
    G3_current_mode = static_cast<int>(byte);
}

void magPro::store_current_direction(char byte) {
    G3_current_direction = static_cast<int>(byte);
}

void magPro::store_waveform(char byte) {
    G3_current_waveform = static_cast<int>(byte);
}

void magPro::store_current_burst_pulses(int pulses) {
    G3_current_burst_pulses = pulses;
}

void magPro::store_current_ipi_index(const std::string& ipi_index) {
    if (ipi_index.size() == 2) {
        int ipi_value = (static_cast<unsigned char>(ipi_index[0]) << 8) | static_cast<unsigned char>(ipi_index[1]);
        G3_current_ipi = static_cast<float>(ipi_value) / 10.0f; // Assuming the value is in 0.1 ms increments
    } else {
        std::cerr << "Invalid IPI index size" << std::endl;
    }
}

void magPro::store_current_ba_ratio(char byte) {
    G3_current_ba_ratio = static_cast<float>(byte) / 10.0f; // Assuming the value is in 0.1 increments
}

void magPro::store_current_enabled(char byte) {
    G3_current_enabled = static_cast<int>(byte);
}

void magPro::store_model(char byte) {
    G3_model = static_cast<int>(byte);
}

void magPro::store_max_pps(int pps) {
    G3_max_pps = pps;
}

void magPro::store_supports_biphasic_burst(char byte) {
    G3_supports_biphasic_burst = static_cast<int>(byte);
}

void magPro::store_mep_values(const std::string& data) {
    if (data.size() == 13 || data.size() == 17) {
        G3_mep_max_amp_value_uv = (static_cast<unsigned char>(data[1]) << 8) | static_cast<unsigned char>(data[2]);
        G3_mep_min_amp_value_uv = (static_cast<unsigned char>(data[3]) << 8) | static_cast<unsigned char>(data[4]);
        G3_mep_max_time_value_us = (static_cast<unsigned char>(data[5]) << 8) | static_cast<unsigned char>(data[6]);
        G3_mep_peak_to_peak_amp_value_uv = (static_cast<unsigned char>(data[7]) << 8) | static_cast<unsigned char>(data[8]);
        G3_mep_timestamp_ms = get_epochtime_ms();
    } else {
        std::cerr << "Invalid MEP data size" << std::endl;
    }
}

void magPro::print_formatted_data_if_all_available() {
    std::cout << "Amplitude A: " << G3_A_amp_percent << "%" << std::endl;
    std::cout << "Amplitude B: " << G3_B_amp_percent << "%" << std::endl;
    std::cout << "di/dt A: " << G3_A_di_dt << std::endl;
    std::cout << "di/dt B: " << G3_B_di_dt << std::endl;
    std::cout << "MEP max amplitude: " << G3_mep_max_amp_value_uv << " uV" << std::endl;
    std::cout << "MEP min amplitude: " << G3_mep_min_amp_value_uv << " uV" << std::endl;
    std::cout << "MEP max time: " << G3_mep_max_time_value_us << " us" << std::endl;
    std::cout << "MEP peak-to-peak amplitude: " << G3_mep_peak_to_peak_amp_value_uv << " uV" << std::endl;
    std::cout << "MEP timestamp: " << G3_mep_timestamp_ms << " ms" << std::endl;
}

void magPro::read_from_serial(unsigned char* buffer, std::size_t size) {
    boost::asio::read(serial, boost::asio::buffer(buffer, size));
}







// Byte string functions
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