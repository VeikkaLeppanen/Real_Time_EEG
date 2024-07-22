#ifndef MAGPRO_H
#define MAGPRO_H

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <boost/asio.hpp>
#include <sys/ioctl.h>



class magPro {
public:
    magPro()
        : serial(io)
    { 
        latest_trigger_time = std::chrono::system_clock::now();
    }

    int connectTriggerPort();
    void trig();
    
    int wait_for_next_package_from_G3();
    void print_formatted_data_if_all_available();

    void set_enable(bool status);
    void set_amplitude(int amplitude);
    void setTriggerTimeLimit(double value) { time_limit = std::max(min_time_limit, std::min(max_time_limit, value)); }
    double getTriggerTimeLimit() { return time_limit; }

    void set_mode(int mode = 0, int direction = 0, int waveform = 1, int burst_pulses = 5, int ipi = 1, double ba_ratio = 1.0, bool delay = true);
    void handle_input_queue();
    bool package_available();
    
    enum class CmdType {
        STATUS = 0,
        AMPLITUDE = 1,
        DI_DT = 2,
        TEMPERATURE = 3,
        MEP_MIN_MAX = 4,
        STATUS_2 = 5,
        ORIG_AMPLITUDE = 6,
        AMPLITUDE_FACTOR = 7,
        PAGE_TRAIN_RUNNING_STATUS = 8,
        MODE = 9
    };

    std::string cmd_type_name(CmdType type) {
        switch (type) {
            case CmdType::STATUS: return "STATUS";
            case CmdType::AMPLITUDE: return "AMPLITUDE";
            case CmdType::DI_DT: return "DI_DT";
            case CmdType::TEMPERATURE: return "TEMPERATURE";
            case CmdType::MEP_MIN_MAX: return "MEP_MIN_MAX";
            case CmdType::STATUS_2: return "STATUS_2";
            case CmdType::ORIG_AMPLITUDE: return "ORIG_AMPLITUDE";
            case CmdType::AMPLITUDE_FACTOR: return "AMPLITUDE_FACTOR";
            case CmdType::PAGE_TRAIN_RUNNING_STATUS: return "PAGE_TRAIN_RUNNING_STATUS";
            case CmdType::MODE: return "MODE";
            default: return "UNKNOWN";
        }
    }
    
private:
    long long get_epochtime_ms();
    std::string get_timestamp_str();
    void cmd_length_4__evaluate_byte_6(char byte);
    void cmd_length_10__mode_9__evaluate_byte_5(char byte);
    void store_mode(char byte);
    void store_current_direction(char byte);
    void store_waveform(char byte);
    void store_current_burst_pulses(int pulses);
    void store_current_ipi_index(const std::string& ipi_index);
    void store_current_ba_ratio(char byte);
    void store_current_enabled(char byte);
    void store_model(char byte);
    void store_max_pps(int pps);
    void store_supports_biphasic_burst(char byte);
    void store_mep_values(const std::string& data);
    void request_G3_to_send_mode_info();
    void verify_mode(int mode, int direction, int waveform, int burst_pulses, int ipi, double ba_ratio);
    void read_from_serial(unsigned char* buffer, std::size_t size);
    int handle_cmd_length_4(const std::string& received_cmd);
    int handle_cmd_length_10(const std::string& received_cmd);
    int handle_cmd_length_mep(const std::string& received_cmd, int length);

    std::vector<uint8_t> create_trig_cmd_byte_str();
    std::vector<uint8_t> create_enable_cmd_byte_str(bool enable);
    std::vector<uint8_t> create_amplitude_cmd_byte_str(uint8_t amplitudeA_value, uint8_t amplitudeB_value);
    std::vector<uint8_t> create_set_mode_cmd_byte_str(int mode, int direction, int waveform, int burst_pulses, int ipi_x10, int ba_ratio_x100);
    std::vector<uint8_t> create_get_mode_cmd_byte_str();
    uint8_t crc8(const std::vector<uint8_t>& data);

    void sleep(double time) { std::this_thread::sleep_for(std::chrono::duration<double>(time)); }
    void print_debug(std::string msg) { if(enable_debug) std::cout << msg << std::endl; }


    // Params
    boost::asio::io_service io;
    boost::asio::serial_port serial;

    // Time limit in milliseconds
    double time_limit;
    const double min_time_limit = 0.1;
    const double max_time_limit = 10000;
    std::chrono::time_point<std::chrono::system_clock> latest_trigger_time;

    bool enable_debug = true;

    double sleep_time_trig = 0.1;
    double sleep_time_set_amp = 1.0;

    double sleep_time_set_mode = 12;
    double sleep_time_set_mode_only_change_IPI_and_amplitude = 0.5;

    double sleep_time_set_enable = sleep_time_set_amp;
    double sleep_time = 0.5;
    
    float G3_current_ba_ratio;              // MagPro G3 X100 limits: 0.2 to 5.0
    int G3_current_direction;               // 0 = Normal, 1 = Reverse
    int G3_current_enabled;                 // 1 = Enabled, 0 = Disabled
    float G3_current_ipi;                   // IPI from 1.0ms to 3000ms
    int G3_current_mode;                    // 0 = Standard, 1 = Power, 2 = Twin, 3 = Dual
    int G3_current_burst_pulses;            // 5, 4, 3 or 2
    int G3_current_waveform;                // 0 = Monophasic, 1 = Biphasic, 2 = Half sine, 3 = Biphasic burst
    int G3_current_page;                    // 8 = MEP page
    int G3_train_running;                   // 0 = Not running, 1 = Running

    int G3_model;                           // 0 = R30, 1 = X100, 2 = R30+Option, 3 = X100+Option, 4 = R30+Option+Mono, 5 = MST
    int G3_max_pps;                         // CHECK VALUE TYPE
    long long G3_supports_biphasic_burst;

    int G3_A_amp_percent;
    int G3_B_amp_percent;
    long long G3_amp_timestamp_ms;

    int G3_A_di_dt;
    int G3_B_di_dt;
    long long G3_di_dt_timestamp_ms;

    int G3_mep_max_amp_value_uv;
    int G3_mep_min_amp_value_uv;
    int G3_mep_max_time_value_us;
    int G3_mep_peak_to_peak_amp_value_uv;
    long long G3_mep_timestamp_ms;
};

#endif // MAGPRO_H