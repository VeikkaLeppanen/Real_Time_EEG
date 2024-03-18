#ifndef DATAHANDLER_H
#define DATAHANDLER_H

#include <mutex>
#include <cstddef> // For size_t
#include <iostream>
#include <chrono>
#include <array>
#include <vector>
#include <cmath>
#include <matio.h>

#include "circularEigenBuffer.h" // Ensure this path is correct
#include "GACorrection.h"

class dataHandler {
public:
    // Default constructor
    dataHandler()
                :   channel_count_(0),
                    sampling_rate_(0),
                    simulation_delivery_rate_(0),
                    GACorr_(GACorrection(0, 0, 0)) {};

    dataHandler(int                  channel_count,
                int                  sampling_rate,
                int                  simulation_delivery_rate)
                :   channel_count_(channel_count),
                    sampling_rate_(sampling_rate),
                    simulation_delivery_rate_(simulation_delivery_rate),
                    GACorr_(GACorrection(channel_count, 25, 100))
    { }

    void reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate);
    void reset_handler(int channel_count, int sampling_rate) { reset_handler(channel_count, sampling_rate, sampling_rate); }

    int simulateData_sin();
    int simulateData_mat();

    void addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger);

    Eigen::MatrixXd getDataInOrder(int downSamplingFactor);
    Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor);
    Eigen::VectorXd getTimeStampsInOrder(int downSamplingFactor);
    Eigen::VectorXd getTriggersInOrder(int downSamplingFactor);

    int get_buffer_capacity() { return buffer_capacity_; }
    int get_channel_count() { return channel_count_; }
    
    void printBufferSize();

    // Eigen::MatrixXd readMatFile(const std::string& fileName);

private:
    std::mutex dataMutex;
    Eigen::MatrixXd sample_buffer_;
    Eigen::VectorXd time_stamp_buffer_;
    Eigen::VectorXd trigger_buffer_;
    size_t current_data_index_ = 0;
    int buffer_length_in_seconds_ = 4;
    int buffer_capacity_;
    int channel_count_;
    int sampling_rate_;
    int simulation_delivery_rate_;

    GACorrection GACorr_;

    // Set this according to the gradient length in samples
    int TA_length = 5000;
    // Set this according to the number of gradients to average over
    int GA_average_length = 25;
    int stimulation_tracker = 100000;
};

#endif // DATAHANDLER_H
