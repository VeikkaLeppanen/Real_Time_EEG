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

    void reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate);
    void reset_handler(int channel_count, int sampling_rate) { reset_handler(channel_count, sampling_rate, sampling_rate); }

    bool isReady() { return (handler_state == WAITING_FOR_STOP); }

    int simulateData_sin();
    int simulateData_mat();

    void addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger, const int &SeqNo);

    int getLatestDataInOrder(Eigen::MatrixXd &output, int number_of_samples);
    Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor);
    Eigen::MatrixXd getMultipleChannelDataInOrder(std::vector<int> channel_indices, int number_of_samples);
    Eigen::MatrixXd getBlockChannelDataInOrder(int first_channel_index, int number_of_channels, int number_of_samples);

    Eigen::VectorXd getTimeStampsInOrder(int downSamplingFactor);
    Eigen::VectorXd getTriggersInOrder(int downSamplingFactor);

    int get_buffer_capacity() { return buffer_capacity_; }
    int get_buffer_length_in_seconds() { return buffer_length_in_seconds_; }
    int get_channel_count() { return channel_count_; }

    void seg_GAcorr_params(int TA_length_input, int GA_average_length_input) { 
        TA_length = TA_length_input; 
        GA_average_length = GA_average_length_input;
    }

    int get_TA_length() { return TA_length; }
    int get_GA_average_length() { return GA_average_length; }
    
    void printBufferSize();

    // Eigen::MatrixXd readMatFile(const std::string& fileName);

private:
    HandlerState handler_state = WAITING_FOR_START;

    std::mutex dataMutex;
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

    GACorrection GACorr_;

    // Set this according to the gradient length in samples
    int TA_length = 5000;
    // Set this according to the number of gradients to average over
    int GA_average_length = 25;
    int stimulation_tracker = 100000;
};

#endif // DATAHANDLER_H
