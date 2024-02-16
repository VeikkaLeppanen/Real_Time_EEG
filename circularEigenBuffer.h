#ifndef CIRCULAREIGENBUFFER_H
#define CIRCULAREIGENBUFFER_H

#include <iostream>
#include <Eigen/Dense>
#include <stdexcept> // For std::out_of_range

class circularEigenBuffer {
public:
    circularEigenBuffer(const int number_of_channels, const int capacity) 
        : data_(Eigen::MatrixXd::Zero(number_of_channels + 1, capacity)), 
          number_of_channels_(number_of_channels), 
          capacity_(capacity), 
          currentIndex_(0) {}
          
    void addSamples(Eigen::VectorXd samples);
    void addSamples(Eigen::MatrixXd samples);
    double getSample(int channel, int index) const;
    Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor = 1) const;
    Eigen::MatrixXd getDataInOrder(int downSamplingFactor = 1) const;
    void print(int channel) const;

    int getCapacity() const { return capacity_; }
    int getCurrentIndex() const { return currentIndex_; }
    int rows() const { return data_.rows(); }
    int cols() const { return data_.cols(); }

private:
    Eigen::MatrixXd data_;
    const int number_of_channels_;
    const int capacity_;
    size_t currentIndex_;
};

#endif // CIRCULAREIGENBUFFER_H
