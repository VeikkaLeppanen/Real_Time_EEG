#include "circularEigenBuffer.h"

// Add a signle sample to each channel
void circularEigenBuffer::addSamples(const Eigen::VectorXd &samples) {
    data_.col(currentIndex_) = samples;
    currentIndex_ = (currentIndex_ + 1) % row_capacity_;
}

// Add multiple samples to each channel
void circularEigenBuffer::addSamples(const Eigen::MatrixXd &samples) {
    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(samples.cols(), static_cast<Eigen::Index>(row_capacity_ - currentIndex_));
    
    // Assign samples that fit to the end
    if (fitToEnd > 0) {
        data_.middleCols(currentIndex_, fitToEnd) = samples.leftCols(fitToEnd);
    }
    
    // Calculate remaining samples that need to wrap around
    int overflow = samples.cols() - fitToEnd;
    
    // If there are any samples that didn't fit, wrap them around to the beginning
    if (overflow > 0) {
        data_.leftCols(overflow) = samples.rightCols(overflow);
    }
    currentIndex_ = (currentIndex_ + samples.cols()) % row_capacity_;
}

// Retrieve the sample at a specific index in the rolling buffer
// Index 0 refers to the oldest sample, and capacity-1 refers to the newest sample.
double circularEigenBuffer::getSample(int channel, int index) const {
    if (index >= row_capacity_) {
        throw std::out_of_range("Index out of bounds");
    }
    // Calculate the actual index based on currentIndex
    int actualIndex = (currentIndex_ + index) % row_capacity_;
    return data_(channel, actualIndex);
}

// Retrieves data from the specified channel in chronological order
Eigen::VectorXd circularEigenBuffer::getChannelDataInOrder(int channel_index, int downSamplingFactor) const {
    if (downSamplingFactor == 0) throw std::invalid_argument("downSamplingFactor must be > 0");

    // Calculate the effective size after downsampling
    size_t downSampledSize = (row_capacity_ + downSamplingFactor - 1) / downSamplingFactor;
    Eigen::VectorXd downSampledData(downSampledSize);

    for (size_t i = 0, j = 0; i < row_capacity_; i += downSamplingFactor, ++j) {
        // Calculate the index in the circular buffer accounting for wrap-around
        size_t index = (currentIndex_ + i) % row_capacity_;
        downSampledData(j) = data_(channel_index, index);
    }

    return downSampledData;
}

// Retrieves data form all channels in chronological order
Eigen::MatrixXd circularEigenBuffer::getDataInOrder(int downSamplingFactor) const {
    if (downSamplingFactor == 0) throw std::invalid_argument("downSamplingFactor must be > 0");

    // Calculate the effective size after downsampling
    size_t downSampledSize = (row_capacity_ + downSamplingFactor - 1) / downSamplingFactor;
    Eigen::MatrixXd downSampledData(data_.rows(), downSampledSize);

    for (size_t i = 0, j = 0; i < row_capacity_; i += downSamplingFactor, ++j) {
        // Calculate the index in the circular buffer accounting for wrap-around
        size_t index = (currentIndex_ + i) % row_capacity_;
        downSampledData.col(j) = data_.col(index);
    }

    return downSampledData;
}

// Print the current state of the buffer for debugging
void circularEigenBuffer::print(int channel) const {
    for (int i = 0; i < row_capacity_; ++i) {
        std::cout << getSample(channel, i) << " ";
    }
    std::cout << "\n";
}
