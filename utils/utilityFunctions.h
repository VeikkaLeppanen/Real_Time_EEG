#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <string>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>

void writeMatrixdToCSV(const std::string& filename, const Eigen::MatrixXd& matrix);
void writeMatrixiToCSV(const std::string& filename, const Eigen::MatrixXi& matrix);
Eigen::MatrixXd readCSV(const std::string &file_path);

Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vec);
Eigen::MatrixXd vectorToColumnMatrixd(const std::vector<double>& vec);
Eigen::MatrixXi vectorToColumnMatrixi(const std::vector<int>& vec);
Eigen::MatrixXd complexVectorToMatrix(const std::vector<std::complex<double>>& complexVec);
Eigen::MatrixXd phaseAngleToMatrix(const std::vector<std::complex<double>>& complexVec);

#endif // UTILITYFUNCTIONS_H