#ifndef PROCESSINGFUNCTIONS_H
#define PROCESSINGFUNCTIONS_H

#include <vector>
#include <string>
#include <Eigen/Dense>

// Function to design a simple FIR low-pass filter using the window method
std::vector<double> designLowPassFilter(int numTaps, double Fs, double Fc) {
    std::vector<double> h(numTaps);
    double r = Fc / (Fs / 2);  // Normalized cutoff frequency

    // Apply the window method (Hann window) and sinc function
    for (int i = 0; i < numTaps; ++i) {
        double M = numTaps - 1;
        double n = i - M / 2;
        if (n == 0.0)
            h[i] = 2 * r;
        else
            h[i] = sin(2 * M_PI * r * n) / (M_PI * n) * (0.5 - 0.5 * cos(2 * M_PI * i / M));

        // Apply the window
        h[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / M);
    }

    return h;
}

// Function to apply a simple FIR filter
Eigen::VectorXd applyFIRFilter(const Eigen::VectorXd& data, const std::vector<double>& b) {
    int n = data.size();
    int nb = b.size();
    Eigen::VectorXd result(n);

    // Zero-padding is not considered here for simplicity and speed
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < nb; ++j) {
            if (i - j >= 0) {
                sum += b[j] * data[i - j];
            }
        }
        result[i] = sum;
    }

    return result;
}

// Function to apply a simple FIR filter to each row of an Eigen matrix
Eigen::MatrixXd applyFIRFilterToMatrix(const Eigen::MatrixXd& dataMatrix, const std::vector<double>& b) {
    int numRows = dataMatrix.rows();
    int numCols = dataMatrix.cols();
    int nb = b.size();
    Eigen::MatrixXd result(numRows, numCols);

    // Process each row with the FIR filter
    for (int i = 0; i < numRows; ++i) {
        Eigen::VectorXd currentRow = dataMatrix.row(i);
        Eigen::VectorXd filteredRow(numCols);

        // Zero-padding is not considered here for simplicity and speed
        for (int j = 0; j < numCols; ++j) {
            double sum = 0.0;
            for (int k = 0; k < nb; ++k) {
                if (j - k >= 0) {
                    sum += b[k] * currentRow[j - k];
                }
            }
            filteredRow[j] = sum;
        }
        result.row(i) = filteredRow;
    }

    return result;
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
Eigen::MatrixXd delayEmbed(const Eigen::MatrixXd& X, int step) {
    int n = X.rows();  // Number of variables in X
    int m = X.cols();  // Length of the signal

    // Total rows in the output matrix considering shifts for each step and the original
    int totalChannels = n * (2 * step + 1);
    Eigen::MatrixXd Y(totalChannels, m); // Initialize matrix for the output

    // Fill the central part of Y with the original data
    int originalStartRow = n * step;
    Y.middleRows(originalStartRow, n) = X;

    // Apply shifts to the left and right
    for (int offset = 1; offset <= step; ++offset) {
        int leftStartRow = originalStartRow - offset * n;
        int rightStartRow = originalStartRow + offset * n;

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

    return Y;
}

void removeBCG(const Eigen::MatrixXd& EEG, Eigen::MatrixXd& CWL, Eigen::MatrixXd& EEG_corrected, int delay) {
    Eigen::MatrixXd expCWL;
    if (delay > 0) {
        // Perform delay embedding if delay is positive
        expCWL = delayEmbed(CWL, (1+2*delay));
    } else {
        expCWL = CWL;
    }

    int num_samples = expCWL.rows();
    Eigen::MatrixXd pinvCWL = expCWL.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Eigen::MatrixXd::Identity(num_samples, num_samples));

    Eigen::MatrixXd EEG_fits = Eigen::MatrixXd::Zero(EEG.rows(), EEG.cols());
    for(int i = 0; i < EEG.rows(); i++) {
        Eigen::MatrixXd Betas = EEG.row(i) * pinvCWL;
        EEG_fits.row(i) = Betas * expCWL;
    }

    EEG_corrected = EEG - EEG_fits;
}


#endif // PROCESSINGFUNCTIONS_H