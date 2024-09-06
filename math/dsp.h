#ifndef DSP_H
#define DSP_H

#include <Eigen/Dense>
#include <complex>

#include <fftw3.h>

Eigen::VectorXd hamming(unsigned int N);
Eigen::VectorXcd spectrum(const Eigen::VectorXd& x, const Eigen::VectorXd& W);
Eigen::MatrixXcd specgram_cx(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl);
Eigen::MatrixXd specgram(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl);
Eigen::VectorXd pwelch(const Eigen::VectorXd& x, unsigned int Nfft = 512, unsigned int Noverl = 256, bool doubleSided = false);
double calculateSNR(const Eigen::VectorXd& data, int overlap, int nfft, double fs, double target_freq, double bandwidth);

double ang_diff(double x, double y);
Eigen::VectorXd ang_diff(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

#endif // DSP_H