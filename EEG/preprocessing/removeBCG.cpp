#include "removeBCG.h"

// Function to perform delay embedding of a signal with incremental shifts and edge value padding
void delayEmbed(const Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int step) {
    int n = X.rows();  // Number of variables in X
    int m = X.cols();  // Length of the signal

    // Fill the central part of Y with the original data
    int originalStartRow = n * step;
    Y.middleRows(originalStartRow, n) = X;

    // Apply shifts to the left and right
    #pragma omp parallel for
    for (int offset = 1; offset <= step; ++offset) {
        int leftStartRow = originalStartRow - offset * n;
        int rightStartRow = originalStartRow + offset * n;

        if (leftStartRow < 0 || rightStartRow + n > Y.rows()) {
            std::cerr << "Error: Invalid input matrix dimensions." << std::endl;
        }

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
}

void removeBCG(const Eigen::MatrixXd& EEG, const Eigen::MatrixXd& expCWL, Eigen::MatrixXd& pinvCWL, Eigen::MatrixXd& EEG_corrected) {
    int num_samples = expCWL.rows();
    
    pinvCWL = expCWL.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Eigen::MatrixXd::Identity(num_samples, num_samples));

    #pragma omp parallel for
    for(int i = 0; i < EEG.rows(); i++) {
        EEG_corrected.row(i) = EEG.row(i) - (EEG.row(i) * pinvCWL) * expCWL;
    }
}