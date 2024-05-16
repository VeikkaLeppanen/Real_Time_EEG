#ifndef PROCESSINGFUNCTIONS_H
#define PROCESSINGFUNCTIONS_H

#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>

// Function declarations
void writeMatrixToCSV(const std::string& filename, const Eigen::MatrixXd& matrix);
void computeButterworthCoefficients(int order, double Fs, double Fc, std::vector<double>& b, std::vector<double>& a);
Eigen::MatrixXd applyIIRFilterToMatrix(const Eigen::MatrixXd& dataMatrix, const std::vector<double>& b, const std::vector<double>& a);
std::vector<double> designLowPassFilter(int numTaps, double Fs, double Fc);
Eigen::VectorXd applyFIRFilter(const Eigen::VectorXd& data, const std::vector<double>& b);
Eigen::MatrixXd applyFIRFilterToMatrix(const Eigen::MatrixXd& dataMatrix, const std::vector<double>& b);
std::vector<double> designBandPassFilter(int numTaps, double Fs, double Fc1, double Fc2);
Eigen::MatrixXd applyBandFIRFilterToMatrix(const Eigen::MatrixXd& dataMatrix, const std::vector<double>& b);
void downsample(const Eigen::MatrixXd& input, Eigen::MatrixXd& output, int factor);
Eigen::MatrixXd delayEmbed(const Eigen::MatrixXd& X, int step);
void removeBCG(const Eigen::MatrixXd& EEG, Eigen::MatrixXd& CWL, Eigen::MatrixXd& EEG_corrected, int delay);

#endif // PROCESSINGFUNCTIONS_H
