#ifndef PHASEESTIMATIONFUNCTIONS_H
#define PHASEESTIMATIONFUNCTIONS_H

#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <tuple>

#include "../../math/dsp.h"

#include <fftw3.h>

void getLSFIRCoeffs_9_13Hz(Eigen::VectorXd& coeffs);

Eigen::VectorXd oddExtension(const Eigen::VectorXd& x, int n);
Eigen::VectorXd applyLSFIRFilter(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs);
Eigen::VectorXd zeroPhaseLSFIR(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs);

Eigen::MatrixXd applyLSFIRFilterMatrix_ret(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs);
void applyLSFIRFilterMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs, Eigen::MatrixXd& output);
void zeroPhaseLSFIRMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs, Eigen::MatrixXd& output);

std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> ls(const Eigen::VectorXd& data, int order, const std::string& norm = "biased", bool allow_singularity = true);
std::vector<double> fitAndPredictAR_LeastSquares(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);
std::vector<double> fitAndPredictAR_Burg(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);
std::vector<double> fitAndPredictAR_YuleWalker(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions);

void levinsonDurbin(const Eigen::VectorXd& r, int order, Eigen::VectorXd& a, double& sigma2, Eigen::VectorXd& k);
Eigen::VectorXd levinsonRecursion(const Eigen::VectorXd &toeplitz, const Eigen::VectorXd &y);
std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> aryule_levinson(const Eigen::VectorXd& data, int order, const std::string& norm = "biased", bool allow_singularity = true);

std::pair<int, double> findTargetPhase(const std::vector<std::complex<double>>& hilbert_signal, Eigen::VectorXd& phaseAngles, int sequence_number, int downsampling_factor, int edge, int prediction_limit, int phase_shift, double stimulation_target);

#endif // PHASEESTIMATIONFUNCTIONS_H
