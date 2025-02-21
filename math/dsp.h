#ifndef DSP_H
#define DSP_H

#include <Eigen/Dense>
#include <complex>
#include <iostream>

#include <fftw3.h>

Eigen::VectorXd computeAutocorrelation(const Eigen::VectorXd& data, int maxLag, const std::string& norm = "biased");
std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> aryule(const Eigen::VectorXd& data, int order, const std::string& norm = "biased", bool allow_singularity = true);

std::vector<std::complex<double>> performFFT(const std::vector<double>& data);
std::vector<std::complex<double>> performFFT(const Eigen::VectorXd& data);
std::vector<std::complex<double>> performIFFT(const std::vector<std::complex<double>>& data);
std::vector<std::complex<double>> hilbertTransform(const std::vector<double>& signal);
std::vector<std::complex<double>> hilbertTransform(const Eigen::VectorXd& signal);

Eigen::VectorXd hamming(unsigned int N);
Eigen::VectorXcd spectrum(const Eigen::VectorXd& x, const Eigen::VectorXd& W);
Eigen::MatrixXcd specgram_cx(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl);
Eigen::MatrixXd specgram(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl);
Eigen::VectorXd pwelch(const Eigen::VectorXd& x, unsigned int Nfft = 512, unsigned int Noverl = 256, bool doubleSided = false);
std::tuple<Eigen::VectorXd, Eigen::VectorXd> computePSD(const Eigen::VectorXd& arParams, double noiseVariance, int nfft, double Fs = 2 * M_PI);
double calculateSNR_max(const Eigen::VectorXd& data, int overlap, int nfft, double fs, double target_freq, double bandwidth, Eigen::VectorXd& Pxx_output);
double calculateSNR_mean(const Eigen::VectorXd& data, int overlap, int nfft, double fs, double target_freq, double bandwidth);

double ang_diff(double x, double y);
Eigen::VectorXd ang_diff(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

#endif // DSP_H