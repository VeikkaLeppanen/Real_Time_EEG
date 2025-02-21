#include "phaseEstimationFunctions.h"

void getLSFIRCoeffs_9_13Hz(Eigen::VectorXd& coeffs) {
    coeffs.resize(71);

    coeffs << -5.76355265e-05, -2.21522410e-03, -4.51913817e-03, -6.92228910e-03,
              -9.37181651e-03, -1.18102996e-02, -1.41771575e-02, -1.64102027e-02,
              -1.84473056e-02, -2.02281241e-02, -2.16958482e-02, -2.27989068e-02,
              -2.34925868e-02, -2.37405117e-02, -2.35159328e-02, -2.28027908e-02,
              -2.15965092e-02, -1.99044896e-02, -1.77462891e-02, -1.51534650e-02,
              -1.21690857e-02, -8.84691437e-03, -5.25028185e-03, -1.45067670e-03,
               2.47391390e-03,  6.44086819e-03,  1.03649495e-02,  1.41604263e-02,
               1.77432351e-02,  2.10331238e-02,  2.39557176e-02,  2.64444474e-02,
               2.84422856e-02,  2.99032408e-02,  3.07935666e-02,  3.10926524e-02,
               3.07935666e-02,  2.99032408e-02,  2.84422856e-02,  2.64444474e-02,
               2.39557176e-02,  2.10331238e-02,  1.77432351e-02,  1.41604263e-02,
               1.03649495e-02,  6.44086819e-03,  2.47391390e-03, -1.45067670e-03,
              -5.25028185e-03, -8.84691437e-03, -1.21690857e-02, -1.51534650e-02,
              -1.77462891e-02, -1.99044896e-02, -2.15965092e-02, -2.28027908e-02,
              -2.35159328e-02, -2.37405117e-02, -2.34925868e-02, -2.27989068e-02,
              -2.16958482e-02, -2.02281241e-02, -1.84473056e-02, -1.64102027e-02,
              -1.41771575e-02, -1.18102996e-02, -9.37181651e-03, -6.92228910e-03,
              -4.51913817e-03, -2.21522410e-03, -5.76355265e-05;
}

Eigen::VectorXd oddExtension(const Eigen::VectorXd& x, int n) {
    int dataSize = x.size();
    if (n < 1) {
        return x;
    }

    Eigen::VectorXd extendedData(dataSize + 2 * n);

    // Left extension
    for (int i = 0; i < n; ++i) {
        int idx = std::min(i + 1, dataSize - 1); // Clamp index within bounds
        extendedData(n - 1 - i) = 2 * x(0) - x(idx);
    }

    // Copy original data
    extendedData.segment(n, dataSize) = x;

    // Right extension
    for (int i = 0; i < n; ++i) {
        int idx = std::max(dataSize - 2 - i, 0); // Clamp index within bounds
        extendedData(dataSize + n + i) = 2 * x(dataSize - 1) - x(idx);
    }

    return extendedData;
}


// Function to apply FIR filter using Eigen
Eigen::VectorXd applyLSFIRFilter(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs) {
    int dataSize = data.size();
    int filterSize = coeffs.size();
    Eigen::VectorXd filteredData = Eigen::VectorXd::Zero(dataSize);

    // Perform convolution
    for (int i = 0; i < dataSize; ++i) {
        for (int j = 0; j < filterSize; ++j) {
            if (i - j >= 0) {
                filteredData(i) += data(i - j) * coeffs(j);
            }
        }
    }

    return filteredData;
}

// Zero-phase filtering equivalent to MATLAB's filtfilt
Eigen::VectorXd zeroPhaseLSFIR(const Eigen::VectorXd& data, const Eigen::VectorXd& coeffs) {
    int extensionSize = 3 * coeffs.size() - 1;
    Eigen::VectorXd paddedData = oddExtension(data, extensionSize);

    // Forward filter pass
    Eigen::VectorXd forwardFiltered = applyLSFIRFilter(paddedData, coeffs);
                 
    // Reverse the data for the backward pass
    Eigen::VectorXd reversedData = forwardFiltered.reverse();
                  
    // Backward filter pass
    Eigen::VectorXd backwardFiltered = applyLSFIRFilter(reversedData, coeffs);
                 
    Eigen::VectorXd orderedData = backwardFiltered.reverse();
    
    // Reverse the result to get the final output
    return orderedData.segment(extensionSize, data.size());
}

Eigen::MatrixXd applyLSFIRFilterMatrix_ret(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs) {
    int numRows = data.rows();
    int numCols = data.cols();
    int filterSize = coeffs.size();
    Eigen::MatrixXd filteredData = Eigen::MatrixXd::Zero(numRows, numCols);

    // Perform convolution for each row
    #pragma omp parallel for
    for (int row = 0; row < numRows; ++row) {
        for (int i = 0; i < numCols; ++i) {
            for (int j = 0; j < filterSize; ++j) {
                if (i - j >= 0) {
                    filteredData(row, i) += data(row, i - j) * coeffs(j);
                }
            }
        }
    }

    return filteredData;
}

void applyLSFIRFilterMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs, Eigen::MatrixXd& filteredData) {
    int numRows = data.rows();
    int numCols = data.cols();
    int filterSize = coeffs.size();
    filteredData.setZero();

    // Perform convolution for each row
    #pragma omp parallel for
    for (int row = 0; row < numRows; ++row) {
        for (int i = 0; i < numCols; ++i) {
            for (int j = 0; j < filterSize; ++j) {
                if (i - j >= 0) {
                    filteredData(row, i) += data(row, i - j) * coeffs(j);
                }
            }
        }
    }
}

void zeroPhaseLSFIRMatrix(const Eigen::MatrixXd& data, const Eigen::VectorXd& coeffs, Eigen::MatrixXd& filteredData) {
    // Forward filter pass
    applyLSFIRFilterMatrix(data, coeffs, filteredData);

    // Reverse the data for the backward pass
    Eigen::MatrixXd reversedData = filteredData;
    #pragma omp parallel for
    for (int row = 0; row < reversedData.rows(); ++row) {
        reversedData.row(row) = reversedData.row(row).reverse();
    }

    // Backward filter pass
    Eigen::MatrixXd backwardFiltered = Eigen::MatrixXd::Zero(data.rows(), data.cols());
    applyLSFIRFilterMatrix(reversedData, coeffs, backwardFiltered);

    // Reverse the result to get the final output
    #pragma omp parallel for
    for (int row = 0; row < backwardFiltered.rows(); ++row) {
        filteredData.row(row) = backwardFiltered.row(row).reverse();
    }
}






















// Function to estimate AR coefficients using least squares method
std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> ls(const Eigen::VectorXd& data, int order, const std::string& norm, bool allow_singularity) {

    Eigen::MatrixXd X(data.size() - order, order);
    Eigen::VectorXd y(data.size() - order);

    for (size_t i = 0; i < data.size() - order; ++i) {
        y(i) = data[i + order];
        for (size_t j = 0; j < order; ++j) {
            X(i, j) = data[i + j];
        }
    }

    // Solve the least squares problem using QR decomposition
    Eigen::VectorXd arParams = X.householderQr().solve(y);
    double sigma2 = (X * arParams - y).squaredNorm() / y.size(); // Calculate sigma^2 as the mean squared error
    Eigen::VectorXd k = Eigen::VectorXd::Zero(order);  // Placeholder, no reflection coefficients computed here

    return std::make_tuple(arParams, sigma2, k);
}

// Function to fit an AR model using the least squares method and predict future values
std::vector<double> fitAndPredictAR_LeastSquares(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions) {
    // Estimate AR coefficients using the Yule-Walker method
    auto [arParams, noiseVariance, reflectionCoeffs] = ls(data, modelOrder);

    // Normalize the data
    double originalMean = data.mean();
    double originalStdDev = std::sqrt((data.array() - originalMean).square().sum() / data.size());
    Eigen::VectorXd normalizedData = (data.array() - originalMean) / originalStdDev;

    // Pre-pend the last 'modelOrder' points of the original data to stabilize predictions
    Eigen::VectorXd prePendedData = Eigen::VectorXd::Zero(modelOrder + numPredictions);
    prePendedData.head(modelOrder) = normalizedData.tail(modelOrder);

    // Predict future values based on the AR model
    std::vector<double> predictions(numPredictions);
    for (size_t i = 0; i < numPredictions; ++i) {
        double predictedValue = 0.0;
        for (size_t j = 0; j < modelOrder; ++j) {
            size_t index = modelOrder + i - 1 - j;
            if (index >= 0 && index < prePendedData.size()) {
                predictedValue += arParams(j) * prePendedData(index);
            }
        }
        prePendedData(modelOrder + i) = predictedValue; // Ensure this index is also valid
        predictions[i] = predictedValue * originalStdDev + originalMean; // Normalize back to original scale

    }

    return predictions;
}

// Function to fit an AR model using the Burg method and predict future values
std::vector<double> fitAndPredictAR_Burg(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions) {
    size_t n = data.size();
    Eigen::VectorXd arParams = Eigen::VectorXd::Zero(modelOrder);  // Initialize AR parameters
    Eigen::VectorXd ef = data;  // Forward errors
    Eigen::VectorXd eb = data;  // Backward errors

    // Iteratively estimate AR parameters using the Burg method
    for (size_t k = 0; k < modelOrder; ++k) {
        // Compute the numerator and the denominator for the reflection coefficient
        double num = 0.0, den = 0.0;
        for (size_t j = 0; j < n - k - 1; ++j) {
            num += ef[j + k + 1] * eb[j];
            den += ef[j + k + 1] * ef[j + k + 1] + eb[j] * eb[j];
        }
        double reflectionCoeff = -2.0 * num / den;

        // Update AR parameters using the reflection coefficient
        arParams[k] = reflectionCoeff;
        for (size_t j = 0; j < k / 2; ++j) {
            double temp = arParams[j];
            arParams[j] += reflectionCoeff * arParams[k - j - 1];
            arParams[k - j - 1] += reflectionCoeff * temp;
        }
        if (k % 2 == 0) {
            arParams[k / 2] += arParams[k / 2] * reflectionCoeff;
        }

        // Update forward and backward prediction errors
        for (size_t j = 0; j < n - k - 1; ++j) {
            double temp = ef[j + k + 1];
            ef[j + k + 1] += reflectionCoeff * eb[j];
            eb[j] += reflectionCoeff * temp;
        }
    }

    // Predict future values based on the AR model
    std::vector<double> predictions(numPredictions, 0.0);
    for (size_t i = 0; i < numPredictions; ++i) {
        double predictedValue = 0;
        for (size_t j = 0; j < modelOrder; ++j) {
            if (n - 1 - j + i >= 0 && n - 1 - j + i < n) {
                predictedValue += arParams[j] * data[n - 1 - j + i];
            } else if (n - 1 - j + i >= n) {
                predictedValue += arParams[j] * predictions[n - 1 - j + i - n];
            }
        }
        predictions[i] = predictedValue;
    }

    return predictions;
}

// Function to fit an AR model using the Yule-Walker method and predict future values
std::vector<double> fitAndPredictAR_YuleWalker(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions) {
    // Estimate AR coefficients using the Yule-Walker method
    auto [arParams, noiseVariance, reflectionCoeffs] = aryule(data, modelOrder);

    // Normalize the data
    double originalMean = data.mean();
    double originalStdDev = std::sqrt((data.array() - originalMean).square().sum() / data.size());
    Eigen::VectorXd normalizedData = (data.array() - originalMean) / originalStdDev;

    // Pre-pend the last 'modelOrder' points of the original data to stabilize predictions
    Eigen::VectorXd prePendedData = Eigen::VectorXd::Zero(modelOrder + numPredictions);
    prePendedData.head(modelOrder) = normalizedData.tail(modelOrder);

    // Predict future values based on the AR model
    std::vector<double> predictions(numPredictions);
    for (size_t i = 0; i < numPredictions; ++i) {
        double predictedValue = 0.0;
        for (size_t j = 0; j < modelOrder; ++j) {
            size_t index = modelOrder + i - 1 - j;
            if (index >= 0 && index < prePendedData.size()) {
                predictedValue += arParams(j) * prePendedData(index);
            }
        }
        prePendedData(modelOrder + i) = predictedValue; // Ensure this index is also valid
        predictions[i] = predictedValue * originalStdDev + originalMean; // Normalize back to original scale

    }

    return predictions;
}










void levinsonDurbin(const Eigen::VectorXd& r, int order, Eigen::VectorXd& a, double& sigma2, Eigen::VectorXd& k) {
    Eigen::VectorXd e(order + 1);
    a = Eigen::VectorXd::Zero(order + 1);
    k = Eigen::VectorXd::Zero(order);
    a(0) = 1.0;
    e(0) = r(0);

    for (int m = 1; m <= order; ++m) {
        double numerator = r(m);
        for (int j = 1; j < m; j++) {  // Ensure not to include the current 'm'
            numerator -= a(j) * r(m - j);
        }
        double km = numerator / e(m - 1);
        k(m - 1) = km;

        Eigen::VectorXd a_new = a;  // Make a copy of the current coefficients
        for (int j = 1; j <= m; j++) {
            a_new(j) = a(j) + km * a(m - j);  // Update based on symmetry
        }
        a = a_new;  // Update the coefficients from the new calculated values

        e(m) = (1 - km * km) * e(m - 1);  // Update the prediction error
    }

    sigma2 = e(order);  // The error variance of the final prediction error
}

Eigen::VectorXd levinsonRecursion(const Eigen::VectorXd &toeplitz, const Eigen::VectorXd &y) {
    const int N = y.size();
    // Correct assertion to match 2 * N + 1 size
    assert(toeplitz.size() == 2 * N + 1); // Ensure toeplitz vector includes zero lag at the center

    Eigen::VectorXd x(N);
    Eigen::VectorXd fn(N), en(N); // Use progressively
    fn.setZero();
    en.setZero();

    int midIndex = N; // Center index of the toeplitz vector
    if (toeplitz[midIndex] == 0) {
        std::cerr << "Error: Central value (zero lag) of toeplitz is zero, inversion impossible.";
        return Eigen::VectorXd(); // Return empty vector if inversion is not possible
    }

    // Initialization based on the zero lag correlation (central value of the toeplitz vector)
    fn[0] = en[0] = 1 / toeplitz[midIndex];
    x[0] = fn[0] * y[0];

    for (int n = 1; n < N; n++) {
        double ef = 0, eb = 0;
        for (int i = 0; i < n; i++) {
            ef += fn[i] * toeplitz[midIndex + n - i];
            eb += en[i] * toeplitz[midIndex - n + i];
        }

        if (1 - ef * eb == 0) {
            std::cerr << "Error: Stability condition failed (1 - ef * eb == 0).";
            return Eigen::VectorXd(); // Return empty vector if stability condition fails
        }

        double c = 1 / (1 - ef * eb); // Inversion for normalization
        for (int i = 0; i <= n; i++) {
            fn[i] *= c;
            en[i] *= c;
        }
        fn[n] = -ef * en[n-1] * c;
        en[n] = -eb * fn[n-1] * c;

        double z = y[n];
        for (int i = 0; i < n; i++) z -= x[i] * toeplitz[midIndex + n - i];
        x[n] = z * en[n];
        for (int i = 0; i <= n; i++) x[i] += en[i] * z;
    }

    return x;
}

std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> aryule_levinson(const Eigen::VectorXd& data, int order, const std::string& norm, bool allow_singularity) {
    Eigen::VectorXd autocorrelations = computeAutocorrelation(data, order, norm);

    int N = order;
    Eigen::VectorXd toeplitz(2 * N + 1);
    // Build symmetric Toeplitz vector from autocorrelations
    for (int i = 0; i <= N; i++) {
        toeplitz[N + i] = toeplitz[N - i] = autocorrelations[i];  // Assuming autocorrelations are symmetric
    }

    Eigen::VectorXd y = autocorrelations.segment(1, N);  // y is r[1] to r[order]
    Eigen::VectorXd arParams = levinsonRecursion(toeplitz, y);

    double sigma2 = 0;  // Compute the final prediction error if necessary
    Eigen::VectorXd k(N);  // Reflection coefficients, not calculated in this scenario

    return std::make_tuple(arParams, sigma2, k);
}

std::pair<int, double> findTargetPhase(
    const std::vector<std::complex<double>>& hilbert_signal,
    Eigen::VectorXd& phaseAngles,
    int sequence_number,
    int downsampling_factor,
    int edge,
    int prediction_limit,
    int phase_shift,
    double stimulation_target)
{
    int target_seqNum = 0;
    double target_phase = 0.0;

    size_t N = hilbert_signal.size();

    // Ensure phaseAngles is appropriately sized
    if (phaseAngles.size() != N)
        phaseAngles.resize(N);

    // Compute the wrapped phase angles using std::arg
    for (size_t i = 0; i < N; ++i) {
        phaseAngles(i) = std::arg(hilbert_signal[i]);
    }

    // Unwrap phase angles to avoid discontinuities
    Eigen::VectorXd phaseAngles_temp(N);
    phaseAngles_temp(0) = phaseAngles(0);
    for (size_t i = 1; i < N; ++i) {
        double delta = phaseAngles(i) - phaseAngles(i - 1);
        // Adjust delta to be within [-π, π]
        delta -= 2 * M_PI * std::round(delta / (2 * M_PI));
        // Cumulatively unwrap the phase
        phaseAngles_temp(i) = phaseAngles_temp(i - 1) + delta;
    }

    // Adjust stimulation_target to be close to phaseAngles_temp[edge]
    double phase_diff = phaseAngles_temp(edge) - stimulation_target;
    double n = std::round(phase_diff / (2 * M_PI));
    double unwrapped_stimulation_target = stimulation_target + n * 2 * M_PI;

    // Compute delta_phase
    Eigen::VectorXd delta_phase = phaseAngles_temp.array() - unwrapped_stimulation_target;

    // Ensure prediction_limit does not exceed N
    size_t limit = std::min(static_cast<size_t>(prediction_limit), N);

    // Find the first crossing point
    for (size_t i = edge + 1; i < limit; ++i) {
        if ((delta_phase(i - 1) < 0 && delta_phase(i) >= 0) ||
            (delta_phase(i - 1) > 0 && delta_phase(i) <= 0)) {

            // Compare absolute differences to find the closer index
            double abs_diff_i_minus_1 = std::abs(delta_phase(i - 1));
            double abs_diff_i = std::abs(delta_phase(i));

            size_t closest_index = (abs_diff_i_minus_1 < abs_diff_i) ? (i - 1) : i;

            // Calculate the target sequence number
            target_seqNum = sequence_number + static_cast<int>((closest_index - edge) * downsampling_factor) + phase_shift;

            // Get the target phase angle from phaseAngles_temp
            target_phase = phaseAngles_temp(closest_index);

            return std::make_pair(target_seqNum, target_phase);
        }
    }

    // If no crossing is found, return 0 and 0.0 for phase
    return std::make_pair(target_seqNum, target_phase);
}

// std::pair<int, double> findTargetPhase(
//     const std::vector<std::complex<double>>& hilbert_signal,
//     Eigen::VectorXd& phaseAngles,
//     int sequence_number,
//     int downsampling_factor,
//     int edge,
//     int prediction_limit,
//     int phase_shift,
//     double stimulation_target)
// {
//     int target_seqNum = 0;
//     double target_phase = 0.0;
//     double adjusted_stimulation_target = stimulation_target;

//     // Compute the phase angles using std::arg
//     for (size_t i = 0; i < hilbert_signal.size(); ++i) {
//         phaseAngles(i) = std::arg(hilbert_signal[i]);
//     }

//     Eigen::VectorXd phaseAngles_temp = phaseAngles;

//     // Unwrap phase angles to avoid discontinuities
//     for (size_t i = 1; i < phaseAngles.size(); ++i) {
//         double delta = phaseAngles_temp(i) - phaseAngles_temp(i - 1);
//         if (delta > M_PI)
//             phaseAngles_temp(i) -= 2 * M_PI;
//         else if (delta < -M_PI)
//             phaseAngles_temp(i) += 2 * M_PI;
//     }

//     // Adjust adjusted_stimulation_target to be in the same range as phaseAngles_temp(0)
//     while (adjusted_stimulation_target - phaseAngles_temp(0) > M_PI)
//         adjusted_stimulation_target -= 2 * M_PI;
//     while (adjusted_stimulation_target - phaseAngles_temp(0) < -M_PI)
//         adjusted_stimulation_target += 2 * M_PI;

//     // Adjust phase angles relative to the stimulation target
//     Eigen::VectorXd delta_phase = phaseAngles_temp.array() - adjusted_stimulation_target;

//     // Find the first crossing point
//     for (size_t i = edge; i < prediction_limit; ++i) {
//         if ((delta_phase(i - 1) < 0 && delta_phase(i) >= 0) ||
//             (delta_phase(i - 1) > 0 && delta_phase(i) <= 0)) {

//             // Compare absolute differences to find the closer index
//             double abs_diff_i_minus_1 = std::abs(delta_phase(i - 1));
//             double abs_diff_i = std::abs(delta_phase(i));

//             int closest_index = (abs_diff_i_minus_1 < abs_diff_i) ? (i - 1) : i;

//             // Calculate the target sequence number
//             target_seqNum = sequence_number + ((closest_index - edge) * downsampling_factor) + phase_shift;

//             // Get the target phase angle from phaseAngles_temp
//             target_phase = phaseAngles_temp(closest_index);

//             return std::make_pair(target_seqNum, target_phase);
//         }
//     }

//     // If no crossing is found, return -1 and 0.0 for phase
//     return std::make_pair(target_seqNum, target_phase);
// }



// std::pair<int, double> findTargetPhase(const std::vector<std::complex<double>>& hilbert_signal,
//                     Eigen::VectorXd& phaseAngles,
//                     int sequence_number,
//                     int downsampling_factor,
//                     int edge,
//                     int phase_shift,
//                     double stimulation_target)
// {
//     int target_seqNum = 0;
//     double target_phase = 0.0;
//     double adjusted_stimulation_target = stimulation_target;

//     // Compute the phase angles using std::arg
//     for (size_t i = 0; i < hilbert_signal.size(); ++i) {
//         phaseAngles(i) = std::arg(hilbert_signal[i]);
//     }
    
//     if (phaseAngles(0) > stimulation_target) {
//         adjusted_stimulation_target = stimulation_target + 2 * M_PI;
//     }
    
//     Eigen::VectorXd phaseAngles_temp = phaseAngles;
    
//     // Unwrap phase angles to avoid discontinuities
//     for (size_t i = 1; i < phaseAngles.size(); ++i) {
//         double delta = phaseAngles_temp(i) - phaseAngles_temp(i - 1);
//         if (delta > M_PI)
//             phaseAngles_temp(i) -= 2 * M_PI;
//         else if (delta < -M_PI)
//             phaseAngles_temp(i) += 2 * M_PI;
//     }
    
//     // Adjust phase angles relative to the stimulation target
//     Eigen::VectorXd delta_phase = phaseAngles_temp.array() - adjusted_stimulation_target;
    
//     // Find the first crossing point
//     for (size_t i = edge; i < delta_phase.size(); ++i) {
//         if ((delta_phase(i - 1) < 0 && delta_phase(i) >= 0) ||
//             (delta_phase(i - 1) > 0 && delta_phase(i) <= 0)) {

//             // Compare absolute differences to find the closer index
//             double abs_diff_i_minus_1 = std::abs(delta_phase(i - 1));
//             double abs_diff_i = std::abs(delta_phase(i));

//             int closest_index = (abs_diff_i_minus_1 < abs_diff_i) ? (i - 1) : i;

//             // std::cout << "Target found at unwrappedPhase " << phaseAngles(closest_index) << std::endl;

//             target_seqNum = sequence_number + ((closest_index - edge) * downsampling_factor) + phase_shift;

//             // Get the target phase angle from phaseAngles_temp
//             target_phase = phaseAngles_temp(closest_index);

//             return std::make_pair(target_seqNum, target_phase);
//         }
//     }
    
//     // If no crossing is found, return -1
//     return std::make_pair(target_seqNum, target_phase);
// }