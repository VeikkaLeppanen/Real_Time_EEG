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
#include <stdexcept>

#include <fftw3.h>

// Function declarations
void writeMatrixdToCSV(const std::string& filename, const Eigen::MatrixXd& matrix);
void writeMatrixiToCSV(const std::string& filename, const Eigen::MatrixXi& matrix);
Eigen::MatrixXd readCSV(const std::string &file_path);

Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vec);
Eigen::MatrixXd vectorToColumnMatrixd(const std::vector<double>& vec);
Eigen::MatrixXi vectorToColumnMatrixi(const std::vector<int>& vec);
Eigen::MatrixXd complexVectorToMatrix(const std::vector<std::complex<double>>& complexVec);
Eigen::MatrixXd phaseAngleToMatrix(const std::vector<std::complex<double>>& complexVec);

void downsample(const Eigen::MatrixXd& input, Eigen::MatrixXd& output, int factor);

void delayEmbed(const Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int step);
void removeBCG(const Eigen::MatrixXd& EEG, const Eigen::MatrixXd& expCWL, Eigen::MatrixXd& pinvCWL, Eigen::MatrixXd& EEG_corrected/*, Eigen::VectorXd& betas*/);

// RT filtering (work in progress)
std::vector<double> createLowPassFilter(int M, double fc, double fs);
std::vector<double> createBandPassFilter(int M, double fc_low, double fc_high, double fs);



void getLSFIRCoeffs_0_80Hz(Eigen::VectorXd& coeffs);
void getLSFIRCoeffs_9_13Hz(Eigen::VectorXd& coeffs);
void designFIR_LS(int numTaps, double f1, double f2, double fs, Eigen::VectorXd& coeffs);

Eigen::VectorXd applyLSFIRFilter(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs);
Eigen::VectorXd zeroPhaseLSFIR(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs);

Eigen::MatrixXd applyLSFIRFilterMatrix_ret(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs);
void applyLSFIRFilterMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs, Eigen::MatrixXd& output);
void zeroPhaseLSFIRMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs, Eigen::MatrixXd& output);

void butterLP(int order, double cutoff, double sampleRate, std::vector<double>& a, std::vector<double>& b);
void getBWCoeffs_8Hz_12Hz(Eigen::VectorXd& a, Eigen::VectorXd& b);
Eigen::VectorXd applyFilter(const Eigen::VectorXd& data, const Eigen::VectorXd& a, const Eigen::VectorXd& b);
Eigen::VectorXd zeroPhaseBW(const Eigen::VectorXd& data, const Eigen::VectorXd& a, const Eigen::VectorXd& b);

std::vector<double> fitAndPredictAR_LeastSquares(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);
std::vector<double> fitAndPredictAR_Burg(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);
std::vector<double> fitAndPredictAR_YuleWalker(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);

Eigen::VectorXd computeAutocorrelation(const Eigen::VectorXd& data, int maxLag);


Eigen::VectorXd computeAutocorrelation_levinson(const Eigen::VectorXd& data, int maxLag, const std::string& norm = "biased");
void levinsonDurbin(const Eigen::VectorXd& r, int order, Eigen::VectorXd& a, double& sigma2, Eigen::VectorXd& k);
Eigen::VectorXd levinsonRecursion(const Eigen::VectorXd &toeplitz, const Eigen::VectorXd &y);
std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> aryule_levinson(const Eigen::VectorXd& data, int order, const std::string& norm = "biased", bool allow_singularity = true);
std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> aryule(const Eigen::VectorXd& data, int order, const std::string& norm = "biased", bool allow_singularity = true);
std::vector<double> fitAndPredictAR_YuleWalker_V2(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);

std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> ls(const Eigen::VectorXd& data, int order, const std::string& norm = "biased", bool allow_singularity = true);
std::vector<double> fitAndPredictAR_LeastSquares_V2(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);


std::vector<std::complex<double>> performFFT(const std::vector<double>& data);
std::vector<std::complex<double>> performFFT(const Eigen::VectorXd& data);
std::vector<std::complex<double>> performIFFT(const std::vector<std::complex<double>>& data);
std::vector<std::complex<double>> hilbertTransform(const std::vector<double>& signal);
std::vector<std::complex<double>> hilbertTransform(const Eigen::VectorXd& signal);

int findTargetPhase(const std::vector<std::complex<double>>& hilbert_signal, Eigen::VectorXd& phaseAngles, int sequence_number, int downsampling_factor, int edge, int phase_shift, double stimulation_target);


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


#endif // PROCESSINGFUNCTIONS_H
