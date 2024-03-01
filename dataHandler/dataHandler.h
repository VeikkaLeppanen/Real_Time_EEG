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

class dataHandler {
public:
    // Default constructor
    dataHandler()
                :   channel_count_(0),
                    sampling_rate_(0),
                    simulation_delivery_rate_(0),
                    sample_packet_size_(0) {};

    dataHandler(int                  channel_count,
                int                  sampling_rate,
                int                  simulation_delivery_rate)
                :   channel_count_(channel_count),
                    sampling_rate_(sampling_rate),
                    simulation_delivery_rate_(simulation_delivery_rate),
                    sample_packet_size_(sampling_rate / simulation_delivery_rate)
    { }

    void reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate);
    void reset_handler(int channel_count, int sampling_rate) { reset_handler(channel_count, channel_count, channel_count); }

    int simulateData();
    
    template<typename Derived>
    void addData(const Eigen::MatrixBase<Derived> &samples);

    Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor);
    Eigen::MatrixXd getDataInOrder(int downSamplingFactor);
    int get_short_buffer_capacity() { return short_buffer_capacity_; }
    int get_short_buffer_num_rows() { return short_buffer_.getNumberOfRows(); }
    
    void printBufferSize();

private:
    std::mutex dataMutex;
    circularEigenBuffer short_buffer_;
    circularEigenBuffer long_buffer_;
    bool useLongBuffer = false;
    int short_buffer_length_in_seconds_ = 2;
    int long_buffer_length_in_seconds_ = 10;
    int short_buffer_capacity_;
    int long_buffer_capacity_;
    int channel_count_;
    int sampling_rate_;
    int simulation_delivery_rate_;
    int sample_packet_size_;
};

#endif // DATAHANDLER_H
