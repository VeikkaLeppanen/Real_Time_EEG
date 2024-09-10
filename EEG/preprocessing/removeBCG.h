#ifndef REMOVEBCG_H
#define REMOVEBCG_H

#include <iostream>
#include <Eigen/Dense>

void delayEmbed(const Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int step);
void removeBCG(const Eigen::MatrixXd& EEG, const Eigen::MatrixXd& expCWL, Eigen::MatrixXd& pinvCWL, Eigen::MatrixXd& EEG_corrected/*, Eigen::VectorXd& betas*/);

#endif // REMOVEBCG_H