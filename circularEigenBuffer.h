#ifndef CIRCULAREIGENBUFFER_H
#define CIRCULAREIGENBUFFER_H

#include <iostream>
#include <Eigen/Dense>
#include <stdexcept> // For std::out_of_range

class circularEigenBuffer {
public:
    circularEigenBuffer(const int number_of_rows, const int row_capacity) 
        : data_(Eigen::MatrixXd::Zero(number_of_rows, row_capacity)), 
          number_of_rows_(number_of_rows), 
          row_capacity_(row_capacity), 
          currentIndex_(0) {}
          
    void addSamples(Eigen::VectorXd &samples);
    void addSamples(Eigen::MatrixXd &samples);
    double getSample(int channel, int index) const;
    Eigen::VectorXd getChannelDataInOrder(int channel_index, int downSamplingFactor = 1) const;
    Eigen::MatrixXd getDataInOrder(int downSamplingFactor = 1) const;
    void print(int channel) const;

    int getNumberOfRows() const { return number_of_rows_; }
    int getCapacity() const { return row_capacity_; }
    int getCurrentIndex() const { return currentIndex_; }

private:
    Eigen::MatrixXd data_;
    const int number_of_rows_;
    const int row_capacity_;
    size_t currentIndex_;
};

#endif // CIRCULAREIGENBUFFER_H
