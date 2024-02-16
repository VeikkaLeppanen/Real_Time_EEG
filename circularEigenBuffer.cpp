#include "circularEigenBuffer.h"

// Add a signle sample to each channel
void circularEigenBuffer::addSamples(Eigen::VectorXd samples) {
    data_.col(currentIndex_) = samples;
    currentIndex_ = (currentIndex_ + 1) % capacity_;
}

// Add multiple samples to each channel
void circularEigenBuffer::addSamples(Eigen::MatrixXd samples) {
    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(samples.cols(), static_cast<Eigen::Index>(capacity_ - currentIndex_));
    
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
}

// Retrieve the sample at a specific index in the rolling buffer
// Index 0 refers to the oldest sample, and capacity-1 refers to the newest sample.
double circularEigenBuffer::getSample(int channel, int index) const {
    if (index >= capacity_) {
        throw std::out_of_range("Index out of bounds");
    }
    // Calculate the actual index based on currentIndex
    int actualIndex = (currentIndex_ + index) % capacity_;
    return data_(channel, actualIndex);
}

Eigen::VectorXd circularEigenBuffer::getChannelDataInOrder(int channel_index, int downSamplingFactor) const {
    // Copy the data from the circular buffer into the ordered vector
    Eigen::VectorXd orderedData(capacity_);

    if (currentIndex_ > 0) {
        int firstSegmentSize = capacity_ - currentIndex_;
        
        // Fill the beginning of orderedData with the "tail" part of the circular buffer
        orderedData.segment(0, firstSegmentSize) = data_.row(channel_index).segment(currentIndex_, firstSegmentSize);
        orderedData.segment(firstSegmentSize, currentIndex_) = data_.row(channel_index).head(currentIndex_);

    } else {
        // If currentIndex_ is 0, the data is already in order.
        orderedData = data_.row(channel_index);
    }

    if (downSamplingFactor > 1) {

        Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<>> downsampledVector(
            orderedData.data(), orderedData.size() / downSamplingFactor, Eigen::InnerStride<>(downSamplingFactor));
            
        return Eigen::VectorXd(downsampledVector);
    }
    return orderedData;
}

Eigen::MatrixXd circularEigenBuffer::getDataInOrder(int downSamplingFactor) const {
    
    // If currentIndex is 0, the data is already in order.
    if (currentIndex_ == 0) { return data_; }

    // Create result matrix with the same size as the original
    Eigen::MatrixXd orderedData(data_.rows(), data_.cols());

    // Determine sizes for the split
    int colsRight = data_.cols() - currentIndex_;
    orderedData.leftCols(colsRight) = data_.rightCols(colsRight);
    orderedData.rightCols(currentIndex_) = data_.leftCols(currentIndex_);

    if (downSamplingFactor > 1) {
        // Calculate the number of columns in the downsampled matrix
        int downsampledCols = (orderedData.cols() + downSamplingFactor - 1) / downSamplingFactor;

        // Use Eigen::Map with Stride to create a downsampled matrix view
        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::Stride<0, Eigen::Dynamic>> downsampledMatrix(
        orderedData.data(), orderedData.rows(), downsampledCols, Eigen::Stride<0, Eigen::Dynamic>(orderedData.outerStride(), downSamplingFactor * orderedData.innerStride()));
        
        return Eigen::MatrixXd(downsampledMatrix);
    }

    return orderedData;
}

// Print the current state of the buffer for debugging
void circularEigenBuffer::print(int channel) const {
    for (int i = 0; i < capacity_; ++i) {
        std::cout << getSample(channel, i) << " ";
    }
    std::cout << "\n";
}
