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
                circularEigenBuffer &short_buffer, 
                circularEigenBuffer &long_buffer,
                int                  channel_count,
                int                  sampling_rate,
                int                  delivery_rate)
                :   dataMutex(mutex),   
                    short_buffer_(short_buffer),
                    long_buffer_(long_buffer),
                    short_buffer_capacity_(short_buffer.getCapacity()),
                    long_buffer_capacity_(long_buffer.getCapacity()),
                    channel_count_(channel_count),
                    sampling_rate_(sampling_rate),
                    delivery_rate_(delivery_rate),
                    sample_packet_size_(sampling_rate / delivery_rate)
    {
        // Parameter checks
        assert(sampling_rate >= delivery_rate && 
               "Error: Sampling rate must be higher than or equal to delivery rate");
        
        assert((long_buffer_capacity_ % short_buffer_capacity_ == 0) && 
                "Error: Long buffer capacity must be an exact multiple of the short buffer capacity");
    }

    int simulateData();
    
    template<typename Derived>
    void addData(const Eigen::MatrixBase<Derived> &samples);
    
    void printBufferSize(int channel);

private:
    std::mutex &dataMutex;
    circularEigenBuffer &short_buffer_;
    circularEigenBuffer &long_buffer_;
    const int short_buffer_capacity_;
    const int long_buffer_capacity_;
    const int channel_count_;
    const int sampling_rate_;
    const int delivery_rate_;
    const int sample_packet_size_;
};

#endif // DATAHANDLER_H
