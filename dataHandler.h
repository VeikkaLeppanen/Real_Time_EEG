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
    dataHandler(circularEigenBuffer &buffer, std::mutex &mutex, int buffer_size, int channel_count)
                                :   dataBuffer(buffer),
                                    dataMutex(mutex),
                                    buffer_size_(buffer_size),
                                    channel_count_(channel_count)
    {}

    int simulateData(int sampling_rate);
    void addData(Eigen::VectorXd samples, double time_stamp);
    void addData(Eigen::MatrixXd samples, double time_stamp);
    void printBufferSize(int channel);

private:
    std::mutex &dataMutex;
    circularEigenBuffer &dataBuffer;
    int buffer_size_;
    int channel_count_;
};

#endif // DATAHANDLER_H
