#ifndef PROCESSINGFUNCTIONS_H
#define PROCESSINGFUNCTIONS_H

#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <omp.h>
#include <numeric>
#include <algorithm>

#include <fftw3.h>
#include <armadillo>

// Function declarations
void writeMatrixToCSV(const std::string& filename, const Eigen::MatrixXd& matrix);
Eigen::MatrixXd readCSV(const std::string &file_path);

Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vec);
Eigen::MatrixXd vectorToColumnMatrix(const std::vector<double>& vec);
Eigen::MatrixXd complexVectorToMatrix(const std::vector<std::complex<double>>& complexVec);
Eigen::MatrixXd phaseAngleToMatrix(const std::vector<std::complex<double>>& complexVec);

std::vector<double> designLowPassFilter(int numTaps, double Fs, double Fc);
std::vector<double> designBandPassFilter(int numTaps, double Fs, double Fc1, double Fc2);

Eigen::MatrixXd applyFIRFilterToMatrix(const Eigen::MatrixXd& dataMatrix, const std::vector<double>& b);

void downsample(const Eigen::MatrixXd& input, Eigen::MatrixXd& output, int factor);

void delayEmbed(const Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int step);
void removeBCG(const Eigen::MatrixXd& EEG, const Eigen::MatrixXd& expCWL, Eigen::MatrixXd& pinvCWL, Eigen::MatrixXd& EEG_corrected/*, Eigen::VectorXd& betas*/);

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

void getLSFIRCoeffs_0_80Hz(Eigen::VectorXd& coeffs);
void getLSFIRCoeffs_9_13Hz(Eigen::VectorXd& coeffs);
void designFIR_LS(int numTaps, double f1, double f2, double fs, Eigen::VectorXd& coeffs);

Eigen::VectorXd applyLSFIRFilter(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs);
Eigen::VectorXd zeroPhaseLSFIR(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs);

Eigen::MatrixXd applyLSFIRFilterMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs);
Eigen::MatrixXd zeroPhaseLSFIRMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs);

void butterLP(int order, double cutoff, double sampleRate, std::vector<double>& a, std::vector<double>& b);
void getBWCoeffs_8Hz_12Hz(Eigen::VectorXd& a, Eigen::VectorXd& b);
Eigen::VectorXd applyFilter(const Eigen::VectorXd& data, const Eigen::VectorXd& a, const Eigen::VectorXd& b);
Eigen::VectorXd zeroPhaseBW(const Eigen::VectorXd& data, const Eigen::VectorXd& a, const Eigen::VectorXd& b);

std::vector<double> fitAndPredictAR_LeastSquares(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);
std::vector<double> fitAndPredictAR_Burg(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);
std::vector<double> fitAndPredictAR_YuleWalker(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);

Eigen::VectorXd computeAutocorrelation(const Eigen::VectorXd& data, int maxLag);

std::vector<std::complex<double>> performFFT(const std::vector<double>& data);
std::vector<std::complex<double>> performFFT(const Eigen::VectorXd& data);
std::vector<std::complex<double>> performIFFT(const std::vector<std::complex<double>>& data);
std::vector<std::complex<double>> hilbertTransform(const std::vector<double>& signal);
std::vector<std::complex<double>> hilbertTransform(const Eigen::VectorXd& signal);




#endif // PROCESSINGFUNCTIONS_H
