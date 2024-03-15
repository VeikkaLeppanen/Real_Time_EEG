#ifndef DATAHANDLER_H
#define DATAHANDLER_H

#include <mutex>
#include <cstddef> // For size_t
#include <iostream>
#include <chrono>
#include <array>
#include <vector>
#include <cmath>
#include "circularEigenBuffer.h" // Ensure this path is correct
#include "GACorrection.h"

class dataHandler {
public:
    // Default constructor
    dataHandler()
                :   channel_count_(0),
                    sampling_rate_(0),
                    simulation_delivery_rate_(0),
                    sample_packet_size_(0),
                    GACorr_(GACorrection(0, 0, 0)) {};

    dataHandler(int                  channel_count,
                int                  sampling_rate,
                int                  simulation_delivery_rate)
                :   channel_count_(channel_count),
                    sampling_rate_(sampling_rate),
                    simulation_delivery_rate_(simulation_delivery_rate),
                    sample_packet_size_(sampling_rate / simulation_delivery_rate),
                    GACorr_(GACorrection(channel_count, 25, 100))
    { }

    void reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate);
    void reset_handler(int channel_count, int sampling_rate) { reset_handler(channel_count, channel_count, channel_count); }

    int simulateData();
    
    void addData(const Eigen::VectorXd &samples, const Eigen::VectorXd &time_stamps, const Eigen::VectorXd &triggers);
    void addData(const Eigen::MatrixXd &samples, const Eigen::VectorXd &time_stamps, const Eigen::VectorXd &triggers);

    Eigen::MatrixXd getDataInOrder(int downSamplingFactor);
    Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor);
    Eigen::VectorXd getTimeStampsInOrder(int downSamplingFactor);
    Eigen::VectorXd getTriggersInOrder(int downSamplingFactor);

    // void addDataGACorr(const Eigen::VectorXd &samples);

    // Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor);
    // Eigen::MatrixXd getDataInOrder(int downSamplingFactor);
    int get_buffer_capacity() { return buffer_capacity_; }
    int get_channel_count() { return channel_count_; }
    
    void printBufferSize();

private:
    std::mutex dataMutex;
    Eigen::MatrixXd sample_buffer_;
    Eigen::VectorXd time_stamp_buffer_;
    Eigen::VectorXd trigger_buffer_;
    size_t current_data_index_ = 0;
    int buffer_length_in_seconds_ = 2;
    int buffer_capacity_;
    int channel_count_;
    int sampling_rate_;
    int simulation_delivery_rate_;
    int sample_packet_size_;

    GACorrection GACorr_;

    // Set this according to the gradient length in samples
    int TA_length = 1000;
    // Set this according to the number of gradients to average over
    int GA_average_length = 25;
    int stimulation_tracker = 10000;
};

#endif // DATAHANDLER_H
