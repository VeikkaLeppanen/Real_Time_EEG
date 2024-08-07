#include "processingFunctions.h"

// Function for writing Eigen matrix to CSV file
void writeMatrixdToCSV(const std::string& filename, const Eigen::MatrixXd& matrix) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j);
                if (j + 1 < matrix.cols()) file << ","; // Comma for next column
            }
            file << "\n"; // Newline for next row
        }
        file.close();
        std::cout << "Data written to " << filename << " successfully." << std::endl;
    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }
}

// Function for writing Eigen matrix to CSV file
void writeMatrixiToCSV(const std::string& filename, const Eigen::MatrixXi& matrix) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j);
                if (j + 1 < matrix.cols()) file << ","; // Comma for next column
            }
            file << "\n"; // Newline for next row
        }
        file.close();
        std::cout << "Data written to " << filename << " successfully." << std::endl;
    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }
}

// Function for reading matrix form a CSV file
Eigen::MatrixXd readCSV(const std::string &file_path) {
    std::ifstream file(file_path);

    // Check if the file stream is valid
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::vector<double> matrix_entries;
    std::string line;
    int num_rows = 0;
    int num_cols = 0;

    // Read each line from the CSV file
    while (getline(file, line)) {
        std::stringstream line_stream(line);
        std::string cell;
        int current_row_cols = 0;

        // Split the line into comma-delimited elements
        while (getline(line_stream, cell, ',')) {
            matrix_entries.push_back(std::stod(cell)); // Convert string to double
            ++current_row_cols;
        }

        if (num_rows == 0) {
            num_cols = current_row_cols; // Set number of columns based on the first row
        }

        ++num_rows;
    }

    file.close();

    // Create an Eigen matrix from the vector of values
    Eigen::MatrixXd matrix(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            matrix(i, j) = matrix_entries[i * num_cols + j];
        }
    }

    return matrix;
}

// Function to convert Eigen::VectorXd to Eigen::MatrixXd for CSV writing
Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vec) {
    // Map the vector as a matrix with one column and vec.size() rows
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), vec.size(), 1);
}

// Function to convert std::vector<double> to a single-column Eigen::MatrixXd
Eigen::MatrixXd vectorToColumnMatrixd(const std::vector<double>& vec) {
    // Map the vector as a single-column matrix
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), vec.size(), 1);
}

// Function to convert std::vector<double> to a single-column Eigen::MatrixXd
Eigen::MatrixXi vectorToColumnMatrixi(const std::vector<int>& vec) {
    // Map the vector as a single-column matrix
    return Eigen::Map<const Eigen::MatrixXi>(vec.data(), vec.size(), 1);
}

// Function to convert std::vector<std::complex<double>> to Eigen::MatrixXd with real and imaginary parts in separate columns
Eigen::MatrixXd complexVectorToMatrix(const std::vector<std::complex<double>>& complexVec) {
    Eigen::MatrixXd combinedMatrix(complexVec.size(), 2);

    for (size_t i = 0; i < complexVec.size(); ++i) {
        combinedMatrix(i, 0) = complexVec[i].real();
        combinedMatrix(i, 1) = complexVec[i].imag();
    }

    return combinedMatrix;
}

// Function to calculate phase angle matrix from a complex vector
Eigen::MatrixXd phaseAngleToMatrix(const std::vector<std::complex<double>>& complexVec) {
    Eigen::MatrixXd phaseAngles(complexVec.size(), 1);

    for (size_t i = 0; i < complexVec.size(); ++i) {
        phaseAngles(i, 0) = std::arg(complexVec[i]);
    }

    return phaseAngles;
}





















// REAL-TIME FILTERING

// Function to create a low-pass windowed-sinc filter kernel
std::vector<double> createLowPassFilter(int M, double fc, double fs) {
    std::vector<double> h(M + 1);
    double sum = 0.0;
    int halfM = M / 2;
    double PI = M_PI;

    for (int i = 0; i <= M; ++i) {
        int n = i - halfM;
        if (n == 0) {
            h[i] = 2 * fc / fs;
        } else {
            h[i] = std::sin(2 * PI * fc / fs * n) / (n * PI);
        }
        // Apply Hamming window
        h[i] *= (0.54 - 0.46 * std::cos(2 * PI * i / M));
        sum += h[i];
    }

    // Normalize the kernel for unity gain at DC
    for (int i = 0; i <= M; ++i) {
        h[i] /= sum;
    }

    return h;
}

// Function to create a bandpass filter using the difference of two low-pass filters
std::vector<double> createBandPassFilter(int M, double fc_low, double fc_high, double fs) {
    auto h_low = createLowPassFilter(M, fc_low, fs);
    auto h_high = createLowPassFilter(M, fc_high, fs);
    std::vector<double> h_band(M + 1);

    for (int i = 0; i <= M; ++i) {
        h_band[i] = h_high[i] - h_low[i];
    }

    return h_band;
}















// OTHER FUNCTIONS

// Downsampling function
void downsample(const Eigen::MatrixXd& input, Eigen::MatrixXd& output, int factor) {
    if (factor <= 0) {
        throw std::invalid_argument("Downsampling factor must be greater than zero.");
    }

    // Compute the number of columns in the downsampled matrix
    int newCols = (input.cols() + factor - 1) / factor;

    for (int j = 0, col = 0; j < input.cols() && col < newCols; j += factor, ++col) {
        output.col(col) = input.col(j);
    }
}

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







// Least squares FIR

// Function to create a Least Squares FIR filter
void designFIR_LS(int numTaps, double f1, double f2, double fs, Eigen::VectorXd& coeffs) {
    int M = numTaps;
    Eigen::MatrixXd A(M, M);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(M);
    double w1 = 2 * M_PI * f1 / fs;  // Lower cutoff frequency (radians)
    double w2 = 2 * M_PI * f2 / fs;  // Upper cutoff frequency (radians)

    // Construct the matrix A and vector b
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            if (i == 0) {
                A(i, j) = w2 / M_PI - w1 / M_PI;
            } else {
                double n = i;
                A(i, j) = (sin(w2 * (j - M/2 + 0.5) * n) - sin(w1 * (j - M/2 + 0.5) * n)) / (M_PI * n);
            }
        }
        b(i) = (i == 0) ? (w2 - w1) : 0;
    }

    // Solve the linear system A * x = b
    coeffs = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    // coeffs = A.colPivHouseholderQr().solve(b);
}


void MultiChannelRealTimeFilter::reset_filter(int numChannels) {
    getLSFIRCoeffs_0_80Hz(filterCoeffs);
    M = filterCoeffs.size();
    buffers = Eigen::MatrixXd::Zero(M, numChannels);
    filteredSamples = Eigen::VectorXd::Zero(numChannels);
}

// Process a new sample vector where each element is the current sample for a channel
Eigen::VectorXd MultiChannelRealTimeFilter::processSample(const Eigen::VectorXd& newSamples) {
    for (int ch = 0; ch < newSamples.size(); ++ch) {
        // Move existing buffer data up one row
        buffers.col(ch).tail(M - 1) = buffers.col(ch).head(M - 1);
        buffers.col(ch)(0) = newSamples(ch);

        // Apply the filter using Eigen dot product
        filteredSamples(ch) = buffers.col(ch).dot(filterCoeffs);
    }
    return filteredSamples;
}

// Function that returns fixed butterworth coefficients for a band pass of 0-80Hz
void getLSFIRCoeffs_0_80Hz(Eigen::VectorXd& coeffs) {
    coeffs.resize(81);

    coeffs << 0.00044545,  0.00066975,  0.00091198,  0.00115233,  0.00136596,  0.00152421,
              0.00159619,  0.00155095,  0.00136009,  0.00100044,  0.00045701, -0.00027433,
             -0.00118445, -0.00224915, -0.00342807, -0.00466456, -0.00588655, -0.00700854,
             -0.00793465, -0.00856266, -0.00878893, -0.00851389, -0.00764786, -0.00611694,
             -0.00386845, -0.00087577,  0.00285786,  0.00729689,  0.01237248,  0.01798326,
              0.02399777,  0.03025859,  0.03658803,  0.04279503,  0.04868311,  0.05405879,
              0.05874021,  0.06256537,  0.06539963,  0.06714197,  0.06772982,  0.06714197,
              0.06539963,  0.06256537,  0.05874021,  0.05405879,  0.04868311,  0.04279503,
              0.03658803,  0.03025859,  0.02399777,  0.01798326,  0.01237248,  0.00729689,
              0.00285786, -0.00087577, -0.00386845, -0.00611694, -0.00764786, -0.00851389,
             -0.00878893, -0.00856266, -0.00793465, -0.00700854, -0.00588655, -0.00466456,
             -0.00342807, -0.00224915, -0.00118445, -0.00027433,  0.00045701,  0.00100044,
              0.00136009,  0.00155095,  0.00159619,  0.00152421,  0.00136596,  0.00115233,
              0.00091198,  0.00066975,  0.00044545;
}

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






// ButterWorth filtering

// Function to calculate the coefficients of a low-pass Butterworth filter
void butterLP(int order, double cutoff, double sampleRate, std::vector<double>& a, std::vector<double>& b) {
    std::vector<std::complex<double>> poles;
    for (int i = 0; i < order; ++i) {
        double theta = M_PI * (2.0 * i + 1.0) / (2.0 * order);
        double realPart = -sin(theta);
        double imagPart = cos(theta);
        std::complex<double> sPole(realPart, imagPart);
        std::complex<double> zPole = (2.0 * sampleRate + sPole) / (2.0 * sampleRate - sPole);  // Bilinear transform
        poles.push_back(zPole);
    }

    std::complex<double> dt = 2.0 * sampleRate;
    std::complex<double> scale = std::pow(std::complex<double>(cutoff * M_PI) / dt, order);

    b.resize(order + 1, 0.0);
    a.resize(order + 1, 0.0);
    b[0] = 1.0;
    a[0] = 1.0;

    for (int i = 0; i < order; ++i) {
        for (int j = i; j >= 0; --j) {
            b[j + 1] += b[j];
            a[j + 1] += real(a[j] * poles[i]);  // Correct use of real() function
        }
    }

    for (int i = 0; i <= order; ++i) {
        a[i] = -real(a[i] * scale);
        b[i] = real(b[i] * scale);
    }

    if (a[0] != 0) {
        for (size_t i = 0; i < a.size(); i++) {
            a[i] /= a[0];
            b[i] /= a[0];
        }
    }
    // Output for debugging
    std::cout << "B (numerator): ";
    for (auto& coeff : b) std::cout << coeff << " ";
    std::cout << "\nA (denominator): ";
    for (auto& coeff : a) std::cout << coeff << " ";
    std::cout << std::endl;
}


// Function that returns fixed butterworth coefficients for a band pass of 8-12Hz
void getBWCoeffs_8Hz_12Hz(Eigen::VectorXd& a, Eigen::VectorXd& b) {
    a.resize(5);
    b.resize(5);

    b << 0.00060985, 0.0, -0.00121971, 0.0, 0.00060985;
    a << 1.0, -3.89919279, 5.73102772, -3.76299529, 0.93138168;
}


// Function to apply the filter in one direction
Eigen::VectorXd applyFilter(const Eigen::VectorXd& data, const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    int n = data.size();
    int m = a.size();
    Eigen::VectorXd filtered(n);

    // Initialize with zero
    filtered.setZero();

    // Direct Form II Transposed implementation
    Eigen::VectorXd z(m - 1);
    z.setZero();

    for (int i = 0; i < n; i++) {
        double out = b(0) * data(i) + z(0);
        for (int j = 1; j < m - 1; j++) {
            z(j - 1) = b(j) * data(i) + z(j) - a(j) * out;
        }
        z(m - 2) = b(m - 1) * data(i) - a(m - 1) * out;
        filtered(i) = out;
    }
    return filtered;
}

// Zero-phase filtering equivalent to MATLAB's filtfilt
Eigen::VectorXd zeroPhaseBW(const Eigen::VectorXd& data, const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    // Forward filter pass
    Eigen::VectorXd forwardFiltered = applyFilter(data, a, b);
    
    // Reverse the data for the backward pass
    Eigen::VectorXd reversedData = forwardFiltered.reverse();
    
    // Backward filter pass
    Eigen::VectorXd backwardFiltered = applyFilter(reversedData, a, b);
    
    // Reverse the result to get the final output
    return backwardFiltered.reverse();
}

















// Function to fit an AR model of given order and predict future values
std::vector<double> fitAndPredictAR_LeastSquares(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions) {
    // Step 1: Prepare the design matrix for the AR model
    Eigen::MatrixXd X(data.size() - modelOrder, modelOrder);
    Eigen::VectorXd y(data.size() - modelOrder);
    
    for (size_t i = 0; i < data.size() - modelOrder; ++i) {
        y(i) = data[i + modelOrder];
        for (size_t j = 0; j < modelOrder; ++j) {
            X(i, j) = data[i + j];
        }
    }

    // Step 2: Estimate AR parameters using linear regression
    Eigen::VectorXd arParams = X.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);

    // Step 3: Predict future values based on the AR model
    std::vector<double> predictions;
    for (size_t i = 0; i < numPredictions; ++i) {
        double predictedValue = 0;
        for (size_t j = 0; j < modelOrder; ++j) {
            size_t index = data.size() - modelOrder + i + j;
            predictedValue += arParams[j] * (index < data.size() ? data(index) : predictions[index - data.size()]);
        }
        predictions.push_back(predictedValue);
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
    // Step 1: Normalize the data
    Eigen::VectorXd normalizedData = data.array() - data.mean();
    normalizedData /= std::sqrt((normalizedData.array() * normalizedData.array()).sum() / normalizedData.size());

    // Step 2: Compute the autocorrelation vector on normalized data
    Eigen::VectorXd autocorrelations = computeAutocorrelation(normalizedData, modelOrder);

    // Step 3: Set up the Yule-Walker equations
    Eigen::MatrixXd R(modelOrder, modelOrder);
    Eigen::VectorXd r(modelOrder);
    for (size_t i = 0; i < modelOrder; ++i) {
        r(i) = autocorrelations(i + 1);
        for (size_t j = 0; j < modelOrder; ++j) {
            R(i, j) = autocorrelations(std::abs(int(i - j)));
        }
    }

    // Solve the Yule-Walker equations to get AR parameters
    Eigen::VectorXd arParams = R.ldlt().solve(r);

    // Step 4: Predict future values based on the AR model
    std::vector<double> predictions(numPredictions);
    for (size_t i = 0; i < numPredictions; ++i) {
        double predictedValue = 0.0;
        for (size_t j = 0; j < modelOrder; ++j) {
            size_t index = data.size() - modelOrder + i + j;
            if (index < data.size()) {
                predictedValue += arParams(j) * normalizedData(index);
            } else {
                predictedValue += arParams(j) * predictions[index - data.size()];
            }
        }
        predictions[i] = predictedValue;
    }

    // Step 5: Rescale the predicted values back to the original data scale
    double originalMean = data.mean();
    double originalStdDev = std::sqrt((data.array() - originalMean).square().sum() / data.size());
    for (size_t i = 0; i < numPredictions; ++i) {
        predictions[i] = predictions[i] * originalStdDev + originalMean;
    }

    return predictions;
}

// Function to compute autocorrelations up to a given lag
Eigen::VectorXd computeAutocorrelation(const Eigen::VectorXd& data, int maxLag) {
    Eigen::VectorXd autocorrelations(maxLag + 1);
    double mean = data.mean();
    double variance = (data.array() - mean).square().sum() / data.size();

    for (int lag = 0; lag <= maxLag; ++lag) {
        if (lag == 0) {
            autocorrelations(lag) = 1;  // Normalized autocorrelation at lag 0
        } else {
            double autocorrelation = 0.0;
            for (int i = 0; i < data.size() - lag; ++i) {
                autocorrelation += (data(i) - mean) * (data(i + lag) - mean);
            }
            autocorrelations(lag) = autocorrelation / ((data.size() - lag) * variance);
        }
    }
    
    return autocorrelations;
}









// Function to compute autocorrelations up to a given lag
Eigen::VectorXd computeAutocorrelation_levinson(const Eigen::VectorXd& data, int maxLag, const std::string& norm) {
    Eigen::VectorXd autocorrelations(maxLag + 1);
    double mean = data.mean();
    double variance = (data.array() - mean).square().sum() / data.size();

    for (int lag = 0; lag <= maxLag; ++lag) {
        if (lag == 0) {
            autocorrelations(lag) = variance;  // Autocorrelation at lag 0 is the variance
        } else {
            double autocorrelation = 0.0;
            for (int i = 0; i < data.size() - lag; ++i) {
                autocorrelation += (data(i) - mean) * (data(i + lag) - mean);
            }
            if (norm == "biased") {
                autocorrelations(lag) = autocorrelation / data.size();
            } else if (norm == "unbiased") {
                autocorrelations(lag) = autocorrelation / (data.size() - lag);
            } else {
                throw std::invalid_argument("Norm must be either 'biased' or 'unbiased'");
            }
        }
    }

    return autocorrelations;
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
    Eigen::VectorXd autocorrelations = computeAutocorrelation_levinson(data, order, norm);

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



// Function to estimate AR coefficients using Yule-Walker method
std::tuple<Eigen::VectorXd, double, Eigen::VectorXd> aryule(const Eigen::VectorXd& data, int order, const std::string& norm, bool allow_singularity) {
    Eigen::VectorXd autocorrelations = computeAutocorrelation_levinson(data, order, norm);

    Eigen::MatrixXd R(order, order);
    Eigen::VectorXd r(order);
    for (int i = 0; i < order; ++i) {
        r(i) = autocorrelations(i + 1);
        for (int j = 0; j < order; ++j) {
            R(i, j) = autocorrelations(std::abs(i - j));
        }
    }

    // Solve the Yule-Walker equations using matrix operations
    Eigen::VectorXd arParams = R.ldlt().solve(r);
    double sigma2 = autocorrelations(0) - arParams.dot(r);
    Eigen::VectorXd k = Eigen::VectorXd::Zero(order);  // Placeholder, no reflection coefficients computed here

    return std::make_tuple(arParams, sigma2, k);
}

// Function to fit an AR model using the Yule-Walker method and predict future values
std::vector<double> fitAndPredictAR_YuleWalker_V2(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions) {
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
std::vector<double> fitAndPredictAR_LeastSquares_V2(const Eigen::VectorXd& data, size_t modelOrder, size_t numPredictions) {
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























// Perform FFT using FFTW
std::vector<std::complex<double>> performFFT(const std::vector<double>& data) {
    int N = data.size();
    std::vector<std::complex<double>> out(N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = data[i];
        in[i][1] = 0;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < N; ++i) {
        out[i] = std::complex<double>(out_fftw[i][0], out_fftw[i][1]);
    }

    fftw_free(in);
    fftw_free(out_fftw);

    return out;
}

// Eigen version of FFT
std::vector<std::complex<double>> performFFT(const Eigen::VectorXd& data) {
    int N = data.size();
    std::vector<std::complex<double>> out(N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = data(i);
        in[i][1] = 0;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < N; ++i) {
        out[i] = std::complex<double>(out_fftw[i][0], out_fftw[i][1]);
    }

    fftw_free(in);
    fftw_free(out_fftw);
    for (int i = 0; i < out.size(); ++i) {
        std::cout << "Freq Bin " << i << ": Magnitude = " << std::abs(out[i]) << std::endl;
    }
    return out;
}   

// Perform inverse FFT using FFTW
std::vector<std::complex<double>> performIFFT(const std::vector<std::complex<double>>& data) {
    int N = data.size();
    std::vector<std::complex<double>> out(N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out_fftw, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < N; ++i) {
        out[i] = std::complex<double>(out_fftw[i][0] / N, out_fftw[i][1] / N); // Normalize by N
    }

    fftw_free(in);
    fftw_free(out_fftw);

    return out;
}

// Apply the Hilbert transform
std::vector<std::complex<double>> hilbertTransform(const std::vector<double>& signal) {
    size_t N = signal.size();
    auto fft_signal = performFFT(signal);

    // Apply the Hilbert transform in the frequency domain
    fft_signal[0] = 0; // DC component zeroed
    for (size_t k = 1; k < N / 2; ++k) {
        fft_signal[k] *= 2;
    }

    if (N % 2 == 0) {
        fft_signal[N / 2] = 0; // Nyquist frequency zeroed if needed
    }

    for (size_t k = N / 2 + 1; k < N; ++k) {
        fft_signal[k] = 0; // Zero the negative frequencies
    }

    // Perform inverse FFT
    return performIFFT(fft_signal);
}

// Eigen version of the Hilbert transform
std::vector<std::complex<double>> hilbertTransform(const Eigen::VectorXd& signal) {
    size_t N = signal.size();
    auto fft_signal = performFFT(signal);

    // Apply the Hilbert transform in the frequency domain
    fft_signal[0] = 0; // DC component zeroed
    for (size_t k = 1; k < N / 2; ++k) {
        fft_signal[k] *= 2;
    }

    if (N % 2 == 0) {
        fft_signal[N / 2] = 0; // Nyquist frequency zeroed if needed
    }

    for (size_t k = N / 2 + 1; k < N; ++k) {
        fft_signal[k] = 0; // Zero the negative frequencies
    }

    // Perform inverse FFT
    return performIFFT(fft_signal);
}


Eigen::VectorXd hamming(unsigned int N) {
    Eigen::VectorXd h(N);
    double a0 = 0.54, a1 = 0.46;
    for (unsigned int i = 0; i < N; ++i) {
        h(i) = a0 - a1 * std::cos(2 * M_PI * i / (N - 1));
    }
    return h;
}


Eigen::VectorXcd spectrum(const Eigen::VectorXd& x, const Eigen::VectorXd& W) {
    int N = x.size();
    Eigen::VectorXcd Pxx(N);

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) {
        in[i][0] = x(i) * W(i); // Real part
        in[i][1] = 0;           // Imaginary part
    }

    fftw_execute(p);

    double wc = W.sum();
    for (int i = 0; i < N; ++i) {
        Pxx(i) = std::complex<double>(out[i][0], out[i][1]) / wc;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return Pxx;
}


Eigen::MatrixXcd specgram_cx(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl) {
    int N = x.size();
    int D = Nfft - Noverl;
    int U = std::max(1, static_cast<int>((N - Nfft) / D + 1)); // Calculate number of columns
    Eigen::MatrixXcd Pw(Nfft, U);
    Eigen::VectorXd W = hamming(Nfft);

    for (int k = 0, m = 0; k <= N - Nfft; k += D, ++m) {
        Eigen::VectorXd xk = x.segment(k, Nfft);
        Pw.col(m) = spectrum(xk, W);
    }

    if (N <= Nfft) {
        Eigen::VectorXd W = hamming(N);
        Pw.resize(N, 1);
        Pw.col(0) = spectrum(x, W);
    }

    return Pw;
}

Eigen::MatrixXd specgram(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl) {
    Eigen::MatrixXcd Pw = specgram_cx(x, Nfft, Noverl);
    return (Pw.array() * Pw.conjugate().array()).real();
}

Eigen::VectorXd pwelch(const Eigen::VectorXd& x, unsigned int Nfft, unsigned int Noverl, bool doubleSided) {
    Eigen::MatrixXd Pxx = specgram(x, Nfft, Noverl);

    if (doubleSided) {
        return Pxx.rowwise().mean();
    } else {
        return Pxx.rowwise().mean().segment(0, Pxx.rows() / 2 + 1);
    }
}

double calculateSNR(const Eigen::VectorXd& data, int overlap, int nfft, double fs, double target_freq, double bandwidth) {
    Eigen::VectorXd psd = pwelch(data.array(), nfft, overlap);

    double df = fs / nfft; // Frequency resolution
    int target_index = static_cast<int>(target_freq / df);
    int half_bandwidth = static_cast<int>(bandwidth / (2 * df));

    double signal_power = 0.0;
    for (int i = target_index - half_bandwidth; i <= target_index + half_bandwidth; ++i) {
        if (i >= 0 && i < psd.size()) {
            signal_power += psd(i);
        }
    }

    double noise_power = 0.0;
    for (int i = 0; i < psd.size(); ++i) {
        if (i < target_index - half_bandwidth || i > target_index + half_bandwidth) {
            noise_power += psd(i);
        }
    }
    noise_power /= (psd.size() - 2 * half_bandwidth);

    return 10 * std::log10(signal_power / noise_power); // SNR in dB
}




int findTargetPhase(const std::vector<std::complex<double>>& hilbert_signal, 
                          Eigen::VectorXd& phaseAngles, 
                          int sequence_number,
                          int downsampling_factor, 
                          int edge, 
                          int phase_shift,
                          double stimulation_target) 
{

    for (std::size_t i = edge; i < hilbert_signal.size(); ++i) {
        phaseAngles(i) = std::arg(hilbert_signal[i]);
        if (i > edge && phaseAngles(i) >= stimulation_target && phaseAngles(i - 1) < stimulation_target) {
            int best_index = std::abs(phaseAngles(i) - stimulation_target) < std::abs(phaseAngles(i - 1) - stimulation_target) ? i : i - 1;
            // std::cout << "Target phase found: " << phaseAngles(best_index) << '\n';
            int trigger_seqNum = sequence_number + (best_index - edge) * downsampling_factor + phase_shift;
            return trigger_seqNum;
        }
    }

    return 0;
}

double ang_diff(double x, double y) {
    std::complex<double> result = std::exp(std::complex<double>(0, x)) / std::exp(std::complex<double>(0, y));
    return std::arg(result);
}

Eigen::VectorXd ang_diff(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    Eigen::VectorXd result(x.size());
    for (int i = 0; i < x.size(); ++i) {
        std::complex<double> num = std::exp(std::complex<double>(0, x[i]));
        std::complex<double> den = std::exp(std::complex<double>(0, y[i]));
        std::complex<double> quotient = num / den;
        result[i] = std::arg(quotient);
    }
    return result;
}