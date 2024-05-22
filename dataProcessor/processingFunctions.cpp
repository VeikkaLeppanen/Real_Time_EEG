#include "processingFunctions.h"

// Function to write Eigen matrix to CSV file
void writeMatrixToCSV(const std::string& filename, const Eigen::MatrixXd& matrix) {
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
    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }
}

// Function to convert Eigen::VectorXd to Eigen::MatrixXd for CSV writing
Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vec) {
    // Map the vector as a matrix with one column and vec.size() rows
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), vec.size(), 1);
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












// FILTERING FUNCTIONS

// Function to design a simple FIR low-pass filter using the window method
std::vector<double> designLowPassFilter(int numTaps, double Fs, double Fc) {
    std::vector<double> h(numTaps);
    double r = Fc / (Fs / 2);  // Normalized cutoff frequency

    // Apply the window method (Hann window) and sinc function
    for (int i = 0; i < numTaps; ++i) {
        double M = numTaps - 1;
        double n = i - M / 2;
        if (n == 0.0)
            h[i] = 2 * r;
        else
            h[i] = sin(2 * M_PI * r * n) / (M_PI * n) * (0.5 - 0.5 * cos(2 * M_PI * i / M));

        // Apply the window
        h[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / M);
    }

    return h;
}

// Function to design a bandpass FIR filter using the window method
std::vector<double> designBandPassFilter(int numTaps, double Fs, double Fc1, double Fc2) {
    std::vector<double> h(numTaps);
    double r1 = Fc1 / (Fs / 2); // Normalized lower cutoff frequency
    double r2 = Fc2 / (Fs / 2); // Normalized upper cutoff frequency

    // Apply the window method (Hann window) and sinc function for bandpass
    for (int i = 0; i < numTaps; ++i) {
        double M = numTaps - 1;
        double n = i - M / 2;
        if (n == 0.0)
            h[i] = 2 * (r2 - r1);
        else
            h[i] = (sin(2 * M_PI * r2 * n) / (M_PI * n)) - (sin(2 * M_PI * r1 * n) / (M_PI * n));

        // Apply the window
        h[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / M);
    }

    return h;
}

// Function to apply a simple FIR filter to each row of an Eigen matrix
Eigen::MatrixXd applyFIRFilterToMatrix(const Eigen::MatrixXd& dataMatrix, const std::vector<double>& b) {
    int numRows = dataMatrix.rows();
    int numCols = dataMatrix.cols();
    int nb = b.size();
    Eigen::MatrixXd result(numRows, numCols);

    // Process each row with the FIR filter
    for (int i = 0; i < numRows; ++i) {
        Eigen::VectorXd currentRow = dataMatrix.row(i);
        Eigen::VectorXd filteredRow(numCols);

        // Zero-padding is not considered here for simplicity and speed
        for (int j = 0; j < numCols; ++j) {
            double sum = 0.0;
            for (int k = 0; k < nb; ++k) {
                if (j - k >= 0) {
                    sum += b[k] * currentRow[j - k];
                }
            }
            filteredRow[j] = sum;
        }
        result.row(i) = filteredRow;
    }

    return result;
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
Eigen::MatrixXd delayEmbed(const Eigen::MatrixXd& X, int step) {
    int n = X.rows();  // Number of variables in X
    int m = X.cols();  // Length of the signal

    // Total rows in the output matrix considering shifts for each step and the original
    int totalChannels = n * (2 * step + 1);
    Eigen::MatrixXd Y(totalChannels, m); // Initialize matrix for the output

    // Fill the central part of Y with the original data
    int originalStartRow = n * step;
    Y.middleRows(originalStartRow, n) = X;

    // Apply shifts to the left and right
    for (int offset = 1; offset <= step; ++offset) {
        int leftStartRow = originalStartRow - offset * n;
        int rightStartRow = originalStartRow + offset * n;

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

    return Y;
}

void removeBCG(const Eigen::MatrixXd& EEG, const Eigen::MatrixXd& CWL, Eigen::MatrixXd& EEG_corrected, int delay) {
    Eigen::MatrixXd expCWL;
    if (delay > 0) {
        // Perform delay embedding if delay is positive
        expCWL = delayEmbed(CWL, (1+2*delay));
    } else {
        expCWL = CWL;
    }

    int num_samples = expCWL.rows();
    Eigen::MatrixXd pinvCWL = expCWL.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Eigen::MatrixXd::Identity(num_samples, num_samples));

    Eigen::MatrixXd EEG_fits = Eigen::MatrixXd::Zero(EEG.rows(), EEG.cols());
    for(int i = 0; i < EEG.rows(); i++) {
        Eigen::MatrixXd Betas = EEG.row(i) * pinvCWL;
        EEG_fits.row(i) = Betas * expCWL;
    }

    EEG_corrected = EEG - EEG_fits;
}












// Perform FFT using FFTW
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