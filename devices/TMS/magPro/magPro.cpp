#include "magPro.h"

int magPro::connectTriggerPort() {

    // Check if the port exists
    // REMEMBER TO GIVE PERMISSIONS TO THE PORT sudo chmod 666 /dev/ttyUSB0
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

    // sleep(sleep_time_trig);
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

// Function to set the mode
void magPro::set_mode(int mode, int direction, int waveform, int burst_pulses, float ipi, float ba_ratio, bool delay) {
    print_debug("set_mode(): Start of function");

    int ipi_x10 = ipi * 10;
    int ba_ratio_x100 = static_cast<int>(round(ba_ratio * 100));

    auto cmd = create_set_mode_cmd_byte_str(mode, direction, waveform, burst_pulses, ipi_x10, ba_ratio_x100);
    boost::asio::write(serial, boost::asio::buffer(cmd));

    print_debug("Mode has been sent.");
    if (delay) {
        print_debug("Sleeping " + std::to_string(sleep_time_set_mode) + " seconds");
        std::this_thread::sleep_for(std::chrono::seconds(static_cast<int>(sleep_time_set_mode)));
    }

    handle_input_queue();
    verify_mode(mode, direction, waveform, burst_pulses, ipi, ba_ratio);
    print_debug("set_mode(): End of function");
}





// Data receiving functions
void magPro::handle_input_queue() {
    while (package_available()) {
        auto received_package = wait_for_next_package_from_G3();
        print_debug("Received package: " + std::to_string(received_package));
    }
}

bool magPro::package_available() {
    int bytes_available;
    if (ioctl(serial.native_handle(), FIONREAD, &bytes_available) == -1) {
        print_debug("Error checking available bytes in serial port");
        return false;
    }
    print_debug("number_of_bytes_in_input_buffer=" + std::to_string(bytes_available));
    return bytes_available > 0;
}

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

void magPro::store_current_ipi_index(const std::string& received_ipi_bytearray) {
    if (received_ipi_bytearray.size() != 2) {
        std::cerr << "Invalid IPI index size" << std::endl;
        return;
    }

    print_debug("received_ipi_bytearray=" + std::to_string(static_cast<int>(received_ipi_bytearray[0])) + " " + std::to_string(static_cast<int>(received_ipi_bytearray[1])));

    int current_ipi_index = static_cast<int>(received_ipi_bytearray[0]) | (static_cast<int>(received_ipi_bytearray[1]) << 8);
    print_debug("current_ipi_index=" + std::to_string(current_ipi_index));

    if (current_ipi_index <= 20) {  // 3000 to 1000 step 100
        G3_current_ipi = 3000 - (100 * current_ipi_index);
    } else if (current_ipi_index > 20 && current_ipi_index <= 30) {  // 1000 to 500 step 50
        G3_current_ipi = 1000 - (50 * (current_ipi_index - 20));
    } else if (current_ipi_index > 30 && current_ipi_index <= 70) {  // 500 to 100 step 10
        G3_current_ipi = 500 - (10 * (current_ipi_index - 30));
    } else if (current_ipi_index > 70 && current_ipi_index <= 150) {  // 100 to 20 step 1.0
        G3_current_ipi = 100 - (current_ipi_index - 70);
    } else if (current_ipi_index > 150 && current_ipi_index <= 170) {  // 20 to 10 step 0.5
        G3_current_ipi = static_cast<int>(round(20 - (0.5 * (current_ipi_index - 150))));
    } else if (current_ipi_index > 170 && current_ipi_index <= 260) {  // 10 to 1 step 0.1
        G3_current_ipi = static_cast<int>(round(10 - (0.1 * (current_ipi_index - 170))));
    } else {
        std::cerr << "current_ipi_index=" + std::to_string(current_ipi_index) + " out of expected range 0 to 270." << std::endl;
        return;
    }

    print_debug("G3_current_ipi=" + std::to_string(G3_current_ipi));
    print_debug("Stored G3_current_ipi =" + std::to_string(G3_current_ipi));
}

void magPro::store_current_ba_ratio(char current_ba_ratio_index) {
    print_debug("_store_current_ba_ratio: current_ba_ratio_index=" + std::to_string(static_cast<int>(current_ba_ratio_index)));
    
    float current_ba_ratio = roundf(5.0f - (0.05f * current_ba_ratio_index) * 100.0f) / 100.0f;
    print_debug("_store_current_ba_ratio: current_ba_ratio=" + std::to_string(current_ba_ratio));

    if (G3_current_mode == 2) { // Assuming TWIN mode is represented by 2
        if (0.2 <= current_ba_ratio && current_ba_ratio <= 5.0) {
            G3_current_ba_ratio = current_ba_ratio;
            print_debug("self.G3_current_ba_ratio=" + std::to_string(G3_current_ba_ratio) + " stored.");
        } else {
            throw std::out_of_range("current_ba_ratio=" + std::to_string(current_ba_ratio) + " out of range. Valid range: 0.2 to 5.0");
        }
    } else {  // Other modes: Standard, Power, Dual
        G3_current_ba_ratio = current_ba_ratio;
        print_debug("self.G3_current_ba_ratio=" + std::to_string(G3_current_ba_ratio) + " stored (unused for current mode).");
    }
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

void magPro::request_G3_to_send_mode_info() {
    print_debug("Running request_G3_to_send_mode_info()");
    
    auto cmd = create_get_mode_cmd_byte_str();
    boost::asio::write(serial, boost::asio::buffer(cmd));
    
    print_debug("Done running request_G3_to_send_mode_info()");
}

void magPro::verify_mode(int mode, int direction, int waveform, int burst_pulses, float ipi, float ba_ratio) {
    print_debug("Requesting G3 to send current mode");
    request_G3_to_send_mode_info();
    print_debug("Done requesting G3 to send response");

    while (wait_for_next_package_from_G3() != static_cast<int>(CmdType::MODE)) {
        std::cout << "Waiting for first MODE package" << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    // Mode data has now been received from G3. Comparing with expected values
    print_debug("G3_current_mode=" + std::to_string(G3_current_mode) + ", mode=" + std::to_string(mode));
    if (G3_current_mode != mode) {
        std::cout << "Warning: Mode could not be verified. Current mode=" + std::to_string(G3_current_mode) + ", requested mode=" + std::to_string(mode) << std::endl;
    }

    print_debug("G3_current_direction=" + std::to_string(G3_current_direction) + ", direction=" + std::to_string(direction));
    if (G3_current_direction != direction) {
        std::cout << "Warning: Direction could not be verified. Current direction=" + std::to_string(G3_current_direction) + ", requested direction=" + std::to_string(direction) << std::endl;
    }

    print_debug("G3_current_waveform=" + std::to_string(G3_current_waveform) + ", waveform=" + std::to_string(waveform));
    if (G3_current_waveform != waveform) {
        std::cout << "Warning: Waveform could not be verified. Current waveform=" + std::to_string(G3_current_waveform) + ", requested waveform=" + std::to_string(waveform) << std::endl;
    }

    print_debug("G3_current_burst_pulses=" + std::to_string(G3_current_burst_pulses) + ", burst_pulses=" + std::to_string(burst_pulses));
    if (G3_current_burst_pulses != burst_pulses) {
        std::cout << "Warning: Burst pulses could not be verified (not relevant for current mode?). Current burst pulses=" + std::to_string(G3_current_burst_pulses) + ", requested burst_pulses=" + std::to_string(burst_pulses) << std::endl;
    }

    print_debug("G3_current_ipi=" + std::to_string(G3_current_ipi) + ", ipi=" + std::to_string(ipi));
    if (G3_current_ipi != ipi) {
        std::cout << "Warning: IPI could not be verified (not relevant for current mode?). Current IPI=" + std::to_string(G3_current_ipi) + ", requested IPI=" + std::to_string(ipi) << std::endl;
    }

    print_debug("G3_current_ba_ratio=" + std::to_string(G3_current_ba_ratio) + ", ba_ratio=" + std::to_string(ba_ratio));
    if (G3_current_ba_ratio != ba_ratio) {
        std::cout << "Warning: B/A ratio could not be verified (not relevant for current mode?). Current B/A ratio=" + std::to_string(G3_current_ba_ratio) + ", requested B/A ratio=" + std::to_string(ba_ratio) << std::endl;
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

// Function to create the command byte string
std::vector<uint8_t> magPro::create_set_mode_cmd_byte_str(int mode, int direction, int waveform, int burst_pulses, int ipi_x10, int ba_ratio_x100) {
    std::vector<uint8_t> cmd_array;

    cmd_array.push_back(0x09); // byte 3: Mode command
    cmd_array.push_back(0x01); // byte 4: 1=set
    cmd_array.push_back(0x00); // byte 5: N/A (only used for "get", not used for "set)

    if (mode < 0 || mode > 3) {
        throw std::out_of_range("'mode' out of range 0-3");
    }
    cmd_array.push_back(mode); // byte 6: mode: 0 (default) = Standard, 1=Power, 2=Twin, 3=Dual

    if (direction < 0 || direction > 1) {
        throw std::out_of_range("'direction' out of range 0-1");
    }
    cmd_array.push_back(direction); // byte 7: direction: 0=normal (default) , 1=reverse

    if (waveform < 0 || waveform > 3) {
        throw std::out_of_range("'waveform' out of range 0-3");
    }
    cmd_array.push_back(waveform); // byte 8: Waveform: 0=Monophasic, 1=Biphasic (default), 2=Half Sine, 3=Biphasic Burst

    if (burst_pulses < 2 || burst_pulses > 5) {
        throw std::out_of_range("'burst_pulses' out of range 2-5");
    }
    cmd_array.push_back(5 - burst_pulses); // byte 9: index 0 = 5 pulses, index 1 = 4 pulses, index 2 = 3 pulses, index 3 = 2 pulses

    if (ipi_x10 < 10 || ipi_x10 > 1000) {
        throw std::out_of_range("'ipi_x10' out of range 10-1000");
    }
    cmd_array.push_back(ipi_x10 / 256); // byte 10: MSB of IPI
    cmd_array.push_back(ipi_x10 % 256); // byte 11: LSB of IPI

    if (ba_ratio_x100 < 20 || ba_ratio_x100 > 500) {
        throw std::out_of_range("'ba_ratio_x100' out of range 20-500");
    }
    cmd_array.push_back(ba_ratio_x100 / 256); // byte 12: MSB of B/A ratio
    cmd_array.push_back(ba_ratio_x100 % 256); // byte 13: LSB of B/A ratio

    uint8_t crc = crc8(cmd_array);
    cmd_array.insert(cmd_array.begin(), 0xfe); // Start byte
    cmd_array.insert(cmd_array.begin() + 1, cmd_array.size()); // Length byte
    cmd_array.push_back(crc); // CRC
    cmd_array.push_back(0xff); // End byte

    return cmd_array;
}

std::vector<uint8_t> magPro::create_get_mode_cmd_byte_str() {
    std::vector<uint8_t> cmd_array = {
        0x09,  // byte 3: Mode command
        0x00,  // byte 4: 0=get
        0x00,  // byte 5: N/A
        0x00,  // byte 6: N/A
        0x00,  // byte 7: N/A
        0x00,  // byte 8: N/A
        0x00,  // byte 9: N/A
        0x00,  // byte 10: N/A
        0x00,  // byte 11: N/A
        0x00,  // byte 12: N/A
        0x00   // byte 13: N/A
    };

    uint8_t crc = crc8(cmd_array);
    std::vector<uint8_t> cmd = {0xfe, static_cast<uint8_t>(cmd_array.size())};
    cmd.insert(cmd.end(), cmd_array.begin(), cmd_array.end());
    cmd.push_back(crc);
    cmd.push_back(0xff);

    return cmd;
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