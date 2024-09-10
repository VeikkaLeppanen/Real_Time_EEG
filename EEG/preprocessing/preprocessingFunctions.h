#ifndef PREPROCESSINGFUNCTIONS_H
#define PREPROCESSINGFUNCTIONS_H

#include <iostream>
#include <Eigen/Dense>

// Real-time filter processor class for multiple channels
class MultiChannelRealTimeFilter {
private:
    Eigen::VectorXd filterCoeffs;
    Eigen::MatrixXd buffers;  // Buffer for each channel
    Eigen::VectorXd filteredSamples;
    int M;

public:
    MultiChannelRealTimeFilter() { }

    void reset_filter(int numChannels);

    // Process a new sample vector where each element is the current sample for a channel
    Eigen::VectorXd processSample(const Eigen::VectorXd& newSamples);
};

void getLSFIRCoeffs_0_80Hz(Eigen::VectorXd& coeffs);


void downsample(const Eigen::MatrixXd& input, Eigen::MatrixXd& output, int factor);

#endif // PREPROCESSINGFUNCTIONS_H