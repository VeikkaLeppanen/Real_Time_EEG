#ifndef DATAHANDLER_H
#define DATAHANDLER_H

#include <mutex>
#include <cstddef> // For size_t
#include <iostream>
#include <chrono>
#include <array>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <Eigen/Dense>

#include "GACorrection.h"
#include "magPro.h"
#include "../dataProcessor/dataProcessor.h"
#include "../dataProcessor/processingFunctions.h"
#include <boost/stacktrace.hpp>

enum HandlerState {
  WAITING_FOR_START,
  WAITING_FOR_STOP
};

class dataHandler {
public:
    // Default constructor
    dataHandler()  // dataProcessor &processor
                :   channel_count_(0),
                    sampling_rate_(0),
                    simulation_delivery_rate_(0),
                    GACorr_(GACorrection(0, 0, 0)),
                    magPro_3G()
    { };

    // Reset functions
    void reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate);
    void reset_handler(int channel_count, int sampling_rate) { reset_handler(channel_count, sampling_rate, sampling_rate); }

    bool isReady() { return (handler_state == WAITING_FOR_STOP); }

    int simulateData_sin();

    // Data handling
    void addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger, const int &SeqNo);

    int getLatestSequenceNumber() { return current_sequence_number_; }
    int getLatestDataInOrder(Eigen::MatrixXd &output, int number_of_samples);
    Eigen::MatrixXd returnLatestDataInOrder(int number_of_samples);
    Eigen::MatrixXd getMultipleChannelDataInOrder(std::vector<int> channel_indices, int number_of_samples);
    Eigen::MatrixXd getBlockChannelDataInOrder(int first_channel_index, int number_of_channels, int number_of_samples);

    Eigen::VectorXd getTimeStampsInOrder(int downSamplingFactor);
    Eigen::VectorXd getTriggersInOrder(int downSamplingFactor);

    int get_buffer_capacity() { return buffer_capacity_; }
    int get_buffer_length_in_seconds() { return buffer_length_in_seconds_; }
    int get_channel_count() { return channel_count_; }
    int get_current_data_index() { return current_data_index_; }

    void setTriggerSource(uint16_t source) { triggerSource = source; }
    void setSourceChannels(std::vector<uint16_t> SourceChannels) {
        source_channels_.resize(SourceChannels.size());
        for (size_t i = 0; i < SourceChannels.size(); ++i) {
            source_channels_(static_cast<int>(i)) = static_cast<int>(SourceChannels[i]);
        }
    }
    Eigen::VectorXi getSourceChannels() const { return source_channels_; }
    
    void setChannelNames(std::vector<std::string> channel_names) { channel_names_ = channel_names; }
    std::vector<std::string> getChannelNames() { return channel_names_; }

    // Gradient artifact correction
    void GACorr_off() { Apply_GACorr = false; }
    void GACorr_on() {
        int stimulation_tracker = 10000000;
        Apply_GACorr = true; 
    }

    void reset_GACorr(int TA_length_input, int GA_average_length_input);
    void reset_GACorr_tracker() { 
        stimulation_tracker = 10000000;
        GACorr_.reset_index();
    }
    int get_TA_length() { return TA_length; }
    int get_GA_average_length() { return GA_average_length; }
    
    // Filtering
    void setFilterState(bool state) { Apply_filter = state; }
    bool getFilterState() { return Apply_filter; }

    // Baseline
    void setBaselineState(bool state) { Apply_baseline = state; }
    bool getBaselineState() { return Apply_baseline; }

    // Triggering
    void setTriggerConnectStatus(bool value) { triggerPortState = value; }
    bool getTriggerConnectStatus() { return triggerPortState; }

    bool getTriggerEnableStatus() { return triggerEnableState; }

    void setTriggerTimeLimit(double value) { magPro_3G.setTriggerTimeLimit(value); }
    double getTriggerTimeLimit() { return magPro_3G.getTriggerTimeLimit(); }

    void insertTrigger(int seqNum) {
        std::lock_guard<std::mutex> lock(triggerMutex);
        triggerSet.insert(seqNum);
    }

    bool shouldTrigger(int seqNum) {
        std::lock_guard<std::mutex> lock(triggerMutex);
        return triggerSet.find(seqNum) != triggerSet.end();
    }

    void removeTrigger(int seqNum) {
        std::lock_guard<std::mutex> lock(triggerMutex);
        triggerSet.erase(seqNum);
    }

    int connectTriggerPort();

    void send_trigger();

    void set_enable(bool status);

    void set_amplitude(int amplitude);

    void magPro_set_mode(int mode = 0, int direction = 0, int waveform = 1, int burst_pulses = 5, float ipi = 1, float ba_ratio = 1.0, bool delay = true);

    void magPro_request_mode_info();

    void get_mode_info(int &mode, int &direction, int &waveform, int &burst_pulses, float &ipi, float &ba_ratio, bool &enabled);

private:
    HandlerState handler_state = WAITING_FOR_START;

    std::mutex dataMutex;
    int triggerSource;
    Eigen::VectorXi source_channels_;
    std::vector<std::string> channel_names_;
    Eigen::MatrixXd sample_buffer_;
    Eigen::VectorXd time_stamp_buffer_;
    Eigen::VectorXd trigger_buffer_;
    size_t current_data_index_ = 0;
    int current_sequence_number_ = 0;
    int buffer_length_in_seconds_ = 20;
    int buffer_capacity_;
    int channel_count_;
    int sampling_rate_;
    int simulation_delivery_rate_;

    // dataProcessor &processor_;

    Eigen::VectorXd processing_sample_vector;

    // Gradient artifact correction
    bool Apply_GACorr = false;
    GACorrection GACorr_;
    int TA_length = 10000;
    int GA_average_length = 25;
    int stimulation_tracker = 10000000;

    // Baseline correction
    bool Apply_baseline = false;
    int baseline_index = 0;
    Eigen::VectorXd baseline_average;
    Eigen::MatrixXd baseline_matrix;
    int baseline_length = 2500;

    // Filter
    bool Apply_filter = false;
    std::vector<double> filterCoeffs_;
    MultiChannelRealTimeFilter RTfilter_;

    // Sending triggers
    magPro magPro_3G;
    bool triggerPortState = false;
    bool triggerEnableState = false;

    std::mutex triggerMutex;
    std::unordered_set<int> triggerSet;
};

#endif // DATAHANDLER_H
