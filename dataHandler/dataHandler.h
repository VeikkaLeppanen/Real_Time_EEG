#ifndef DATAHANDLER_H
#define DATAHANDLER_H

#include <mutex>
#include <cstddef> // For size_t
#include <iostream>
#include <chrono>
#include <array>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "GACorrection.h"
#include "../dataProcessor/dataProcessor.h"

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
                    GACorr_(GACorrection(0, 0, 0)) {};
                    // processor_(processor) {};

    // Reset functions
    void reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate);
    void reset_handler(int channel_count, int sampling_rate) { reset_handler(channel_count, sampling_rate, sampling_rate); }

    bool isReady() { return (handler_state == WAITING_FOR_STOP); }

    int simulateData_sin();
    int simulateData_mat();

    void addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger, const int &SeqNo);

    int getLatestDataInOrder(Eigen::MatrixXd &output, int number_of_samples);
    Eigen::MatrixXd returnLatestDataInOrder(int number_of_samples);
    Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor);
    Eigen::MatrixXd getMultipleChannelDataInOrder(std::vector<int> channel_indices, int number_of_samples);
    Eigen::MatrixXd getBlockChannelDataInOrder(int first_channel_index, int number_of_channels, int number_of_samples);

    Eigen::VectorXd getTimeStampsInOrder(int downSamplingFactor);
    Eigen::VectorXd getTriggersInOrder(int downSamplingFactor);

    int get_buffer_capacity() { return buffer_capacity_; }
    int get_buffer_length_in_seconds() { return buffer_length_in_seconds_; }
    int get_channel_count() { return channel_count_; }

    void setTriggerSource(uint16_t source) { triggerSource = source; }
    void setSourceChannels(std::vector<uint16_t> SourceChannels) {
        source_channels_.resize(SourceChannels.size());
        for (size_t i = 0; i < SourceChannels.size(); ++i) {
            source_channels_(static_cast<int>(i)) = static_cast<int>(SourceChannels[i]);
        }
    }
    Eigen::VectorXi getSourceChannels() { return source_channels_; }
    
    void setChannelNames(std::vector<std::string> channel_names) { channel_names_ = channel_names; }
    std::vector<std::string> getChannelNames() { return channel_names_; }

    // Gradient artifact correction functions
    void GACorr_off() { GACorr_running = false; }
    void GACorr_on() {
        int stimulation_tracker = 10000000;
        GACorr_running = true; 
    }

    void reset_GACorr(int TA_length_input, int GA_average_length_input);
    void reset_GACorr_tracker() { 
        stimulation_tracker = 10000000;
        GACorr_.reset_index();
    }
    int get_TA_length() { return TA_length; }
    int get_GA_average_length() { return GA_average_length; }
    
    void printBufferSize();

    // Eigen::MatrixXd readMatFile(const std::string& fileName);

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

    bool GACorr_running = false;
    GACorrection GACorr_;
    int TA_length = 5000;
    int GA_average_length = 25;
    int stimulation_tracker = 10000000;
};

#endif // DATAHANDLER_H
