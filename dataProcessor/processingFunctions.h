#ifndef PROCESSINGFUNCTIONS_H
#define PROCESSINGFUNCTIONS_H

#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <fftw3.h>
#include <complex>
#include <cmath>

// Function declarations
void writeMatrixToCSV(const std::string& filename, const Eigen::MatrixXd& matrix);
Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vec);
Eigen::MatrixXd complexVectorToMatrix(const std::vector<std::complex<double>>& complexVec);
Eigen::MatrixXd phaseAngleToMatrix(const std::vector<std::complex<double>>& complexVec);


std::vector<double> designLowPassFilter(int numTaps, double Fs, double Fc);
std::vector<double> designBandPassFilter(int numTaps, double Fs, double Fc1, double Fc2);

Eigen::MatrixXd applyFIRFilterToMatrix(const Eigen::MatrixXd& dataMatrix, const std::vector<double>& b);

void downsample(const Eigen::MatrixXd& input, Eigen::MatrixXd& output, int factor);

Eigen::MatrixXd delayEmbed(const Eigen::MatrixXd& X, int step);
void removeBCG(const Eigen::MatrixXd& EEG, const Eigen::MatrixXd& CWL, Eigen::MatrixXd& EEG_corrected, int delay);

// RT filtering (work in progress)
std::vector<double> createLowPassFilter(int M, double fc, double fs);
std::vector<double> createBandPassFilter(int M, double fc_low, double fc_high, double fs);

// Real-time filter processor class for multiple channels
class MultiChannelRealTimeFilter {
private:
    std::vector<double> filterCoeffs;
    std::vector<std::vector<double>> buffers;  // Buffer for each channel
    int M;

public:
    MultiChannelRealTimeFilter() { }

    void reset_filter(const std::vector<double>& coeffs, int numChannels) {
        filterCoeffs = coeffs;
        M = coeffs.size();
        buffers.resize(numChannels, std::vector<double>(M, 0.0));
    }

    // Process a new sample vector where each element is the current sample for a channel
    Eigen::VectorXd processSample(const Eigen::VectorXd& newSamples) {
        Eigen::VectorXd filteredSamples(newSamples.size());
        for (int ch = 0; ch < newSamples.size(); ++ch) {
            // Shift buffer contents to make room for the new sample
            std::rotate(buffers[ch].rbegin(), buffers[ch].rbegin() + 1, buffers[ch].rend());
            buffers[ch][0] = newSamples[ch];

            // Apply the filter
            double filteredSample = 0.0;
            for (int i = 0; i < M; ++i) {
                filteredSample += filterCoeffs[i] * buffers[ch][i];
            }
            filteredSamples[ch] = filteredSample;
        }
        return filteredSamples;
    }
};


std::vector<std::complex<double>> performFFT(const Eigen::VectorXd& data);
std::vector<std::complex<double>> performIFFT(const std::vector<std::complex<double>>& data);
std::vector<std::complex<double>> hilbertTransform(const Eigen::VectorXd& signal);




#endif // PROCESSINGFUNCTIONS_H
