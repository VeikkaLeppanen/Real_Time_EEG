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
    dataHandler(std::mutex          &mutex, 
                circularEigenBuffer &shortBuffer, 
                circularEigenBuffer &longBuffer,
                int                  channel_count,
                int                  sampling_rate,
                int                  delivery_rate)
                :   dataMutex(mutex),   
                    shortBuffer_(shortBuffer),
                    longBuffer_(longBuffer),
                    short_buffer_capacity_(shortBuffer.getCapacity()),
                    long_buffer_capacity_(longBuffer.getCapacity()),
                    channel_count_(channel_count),
                    sampling_rate_(sampling_rate),
                    delivery_rate_(delivery_rate),
                    sample_packet_size_(sampling_rate / delivery_rate)
    {}

    int simulateData();
    
    template<typename Derived>
    void addData(const Eigen::MatrixBase<Derived> &samples);
    
    void printBufferSize(int channel);

private:
    std::mutex &dataMutex;
    circularEigenBuffer &shortBuffer_;
    circularEigenBuffer &longBuffer_;
    int short_buffer_capacity_;
    int long_buffer_capacity_;
    int channel_count_;
    int sampling_rate_;
    int delivery_rate_;
    int sample_packet_size_;
};

#endif // DATAHANDLER_H
