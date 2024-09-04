#include "preprocessingFunctions.h"

void MultiChannelRealTimeFilter::reset_filter(int numChannels) {
    getLSFIRCoeffs_0_80Hz(filterCoeffs);
    M = filterCoeffs.size();
    buffers = Eigen::MatrixXd::Zero(M, numChannels);
    filteredSamples = Eigen::VectorXd::Zero(numChannels);
}

// Process a new sample vector where each element is the current sample for a channel
Eigen::VectorXd MultiChannelRealTimeFilter::processSample(const Eigen::VectorXd& newSamples) {
    for (int ch = 0; ch < newSamples.size(); ++ch) {
        // Move existing buffer data up one row
        buffers.col(ch).tail(M - 1) = buffers.col(ch).head(M - 1);
        buffers.col(ch)(0) = newSamples(ch);

        // Apply the filter using Eigen dot product
        filteredSamples(ch) = buffers.col(ch).dot(filterCoeffs);
    }
    return filteredSamples;
}

// Function that returns fixed butterworth coefficients for a band pass of 0-80Hz
void getLSFIRCoeffs_0_80Hz(Eigen::VectorXd& coeffs) {
    coeffs.resize(81);

    coeffs << 0.00044545,  0.00066975,  0.00091198,  0.00115233,  0.00136596,  0.00152421,
              0.00159619,  0.00155095,  0.00136009,  0.00100044,  0.00045701, -0.00027433,
             -0.00118445, -0.00224915, -0.00342807, -0.00466456, -0.00588655, -0.00700854,
             -0.00793465, -0.00856266, -0.00878893, -0.00851389, -0.00764786, -0.00611694,
             -0.00386845, -0.00087577,  0.00285786,  0.00729689,  0.01237248,  0.01798326,
              0.02399777,  0.03025859,  0.03658803,  0.04279503,  0.04868311,  0.05405879,
              0.05874021,  0.06256537,  0.06539963,  0.06714197,  0.06772982,  0.06714197,
              0.06539963,  0.06256537,  0.05874021,  0.05405879,  0.04868311,  0.04279503,
              0.03658803,  0.03025859,  0.02399777,  0.01798326,  0.01237248,  0.00729689,
              0.00285786, -0.00087577, -0.00386845, -0.00611694, -0.00764786, -0.00851389,
             -0.00878893, -0.00856266, -0.00793465, -0.00700854, -0.00588655, -0.00466456,
             -0.00342807, -0.00224915, -0.00118445, -0.00027433,  0.00045701,  0.00100044,
              0.00136009,  0.00155095,  0.00159619,  0.00152421,  0.00136596,  0.00115233,
              0.00091198,  0.00066975,  0.00044545;
}


// Downsampling function
void downsample(const Eigen::MatrixXd& input, Eigen::MatrixXd& output, int factor) {
    if (factor <= 0) {
        throw std::invalid_argument("Downsampling factor must be greater than zero.");
    }

    // Compute the number of columns in the downsampled matrix
    int newCols = (input.cols() + factor - 1) / factor;

    for (int j = 0, col = 0; j < input.cols() && col < newCols; j += factor, ++col) {
        output.col(col) = input.col(j);
    }
}

// Function to perform delay embedding of a signal with incremental shifts and edge value padding
void delayEmbed(const Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int step) {
    int n = X.rows();  // Number of variables in X
    int m = X.cols();  // Length of the signal

    // Fill the central part of Y with the original data
    int originalStartRow = n * step;
    Y.middleRows(originalStartRow, n) = X;

    // Apply shifts to the left and right
    #pragma omp parallel for
    for (int offset = 1; offset <= step; ++offset) {
        int leftStartRow = originalStartRow - offset * n;
        int rightStartRow = originalStartRow + offset * n;

        if (leftStartRow < 0 || rightStartRow + n > Y.rows()) {
            std::cerr << "Error: Invalid input matrix dimensions." << std::endl;
        }

        // Left shifts (data moves right)
        for (int row = 0; row < n; ++row) {
            Y.block(leftStartRow + row, offset, 1, m - offset) = X.block(row, 0, 1, m - offset);
            // Edge value padding on the left
            Y.block(leftStartRow + row, 0, 1, offset).setConstant(X(row, 0));
        }

        // Right shifts (data moves left)
        for (int row = 0; row < n; ++row) {
            Y.block(rightStartRow + row, 0, 1, m - offset) = X.block(row, offset, 1, m - offset);
            // Edge value padding on the right
            Y.block(rightStartRow + row, m - offset, 1, offset).setConstant(X(row, m - 1));
        }
    }
}

void removeBCG(const Eigen::MatrixXd& EEG, const Eigen::MatrixXd& expCWL, Eigen::MatrixXd& pinvCWL, Eigen::MatrixXd& EEG_corrected) {
    int num_samples = expCWL.rows();
    
    pinvCWL = expCWL.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Eigen::MatrixXd::Identity(num_samples, num_samples));

    #pragma omp parallel for
    for(int i = 0; i < EEG.rows(); i++) {
        EEG_corrected.row(i) = EEG.row(i) - (EEG.row(i) * pinvCWL) * expCWL;
    }
}
