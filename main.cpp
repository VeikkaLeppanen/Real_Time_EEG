#include <vector>
#include <thread>
#include <mutex>
#include <omp.h>
#include <string>
#include <csignal>
#include <chrono>
#include <fstream>

#include "dataProcessor/dataProcessor.h"
#include "dataHandler/dataHandler.h"
#include "eeg_bridge/eeg_bridge.h"

#include "UI/mainwindow.h"
#include <QApplication>
#include <QWidget>

// #include "matplotlibcpp.h"
// namespace plt = matplotlibcpp;

// Simulation loop parameters
const uint8_t CHANNEL_COUNT = 12;
const uint32_t SAMPLING_RATE = 5000;
const uint32_t DELIVERY_RATE = 5000;
const uint32_t DOWNSAMPLING_FACTOR = 1;


// In case of bind failed the previous process can be terminated with the following commands on linux
// lsof -i :8080        Find PID
// kill -9 PID          Kill the process

// This is used to terminate the program with Ctrl+C
volatile std::sig_atomic_t signal_received = 0;
void signal_handler(int signal) {
    signal_received = 1;
}




// DATA AQUISITION AND RAW DATA PLOTTING

// Loop for the data collecting and storing
void dataAcquisitionLoop(dataHandler &handler) {

    EegBridge bridge;
    bridge.bind_socket();
    bridge.spin(handler, signal_received);

    std::cout << "Exiting dataAcquisitionLoop" << '\n';
}







// DATA PROCESSING AND PROCESSED DATA PLOTTING

// Filtering
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

// Function to apply a simple FIR filter
Eigen::VectorXd applyFIRFilter(const Eigen::VectorXd& data, const std::vector<double>& b) {
    int n = data.size();
    int nb = b.size();
    Eigen::VectorXd result(n);

    // Zero-padding is not considered here for simplicity and speed
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < nb; ++j) {
            if (i - j >= 0) {
                sum += b[j] * data[i - j];
            }
        }
        result[i] = sum;
    }

    return result;
}


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

void removeBCG(const Eigen::MatrixXd& EEG, Eigen::MatrixXd& CWL, Eigen::MatrixXd& EEG_corrected, int delay) {
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

// Loop for the data processing
void dataProcessingLoop(dataHandler &handler) {
    
    // Initializing parameters and matrices for algorithms
    int n_eeg_channels = 5;
    int n_cwl_channels = 7;
    int n_samples = 10000;
    Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(handler.get_channel_count(), n_samples);
    Eigen::MatrixXd EEG = Eigen::MatrixXd::Zero(n_eeg_channels, n_samples);
    Eigen::MatrixXd CWL = Eigen::MatrixXd::Zero(n_cwl_channels, n_samples);

    // Filtering
    double Fs = 1000;  // Sampling frequency
    double Fc = 120;   // Desired cutoff frequency
    int numTaps = 51;  // Length of the FIR filter
    std::vector<double> filterCoeffs = designLowPassFilter(numTaps, Fs, Fc);
    Eigen::MatrixXd EEG_filtered = Eigen::MatrixXd::Zero(n_eeg_channels, n_samples);
    Eigen::MatrixXd CWL_filtered = Eigen::MatrixXd::Zero(n_cwl_channels, n_samples);

    // Downsampling
    int downsampling_factor = 10;
    int newCols = (n_samples + downsampling_factor - 1) / downsampling_factor;
    Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_eeg_channels, newCols);
    Eigen::MatrixXd CWL_downsampled = Eigen::MatrixXd::Zero(n_cwl_channels, newCols);

    // BCG correction
    int delay = 0;
    Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_eeg_channels, newCols);

    // TESTING MATRICES
    Eigen::VectorXd Betas_temp = Eigen::VectorXd::Zero(n_cwl_channels);
    Eigen::MatrixXd Betas = Eigen::MatrixXd::Zero(127, n_cwl_channels);
    Eigen::MatrixXd EEG_output = Eigen::MatrixXd::Zero(127, newCols);
    double total_time = 0.0;
    int count = 0;


    while (!handler.isReady()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    while (!signal_received) {

        // Copy data from wanted EEG channels and CWL channels
        int sequence_number = handler.getLatestDataInOrder(all_channels, n_samples);
        EEG = all_channels.middleRows(0, n_eeg_channels);
        CWL = all_channels.middleRows(5, n_cwl_channels);

        // if (sequence_number > 0 && sequence_number % 10000 == 0) {
        // std::cout << "SeqNum: " << sequence_number << '\n';

        // LOWPASS filter 120Hz   BANDPASS 0.33-125Hz FIR
        EEG_filtered = applyFIRFilter(EEG, filterCoeffs);
        CWL_filtered = applyFIRFilter(CWL, filterCoeffs);

        // DONWSAMPLING
        downsample(EEG_filtered, EEG_downsampled, downsampling_factor);
        downsample(CWL_filtered, CWL_downsampled, downsampling_factor);

        // auto start = std::chrono::high_resolution_clock::now();

        // removeBCG
        removeBCG(EEG_downsampled, CWL_downsampled, EEG_corrected, delay);

        // std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
        // std::cout << "Time taken: " << elapsed.count() << " seconds." << std::endl;
        // total_time += elapsed.count();
        // count++;
    }

    // if (count > 0) {
    //     double average_time = total_time / count;
    //     std::cout << "Total time taken: " << total_time << " seconds." << std::endl;
    //     std::cout << "Average time taken: " << average_time << " seconds." << std::endl;

    //     writeMatrixToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/Betas.csv", Betas);
    //     writeMatrixToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/EEG_corrected.csv", EEG_output);
    // }

    std::cout << "Exiting dataProcessingLoop" << '\n';
}








// int main(int argc, char *argv[])
// {
//     dataHandler handler;

//     QApplication a(argc, argv);
//     MainWindow w(handler, signal_received);
//     w.show();
//     return a.exec();
// }




int main() {
    std::signal(SIGINT, signal_handler);
    
    // size_t thread_count = 6;
    // dataProcessor processor(thread_count);
    // dataHandler handler(processor);
    dataHandler handler;

    // std::thread dataThread(dataSimulationLoop, std::ref(handler));
    std::thread dataThread(dataAcquisitionLoop, std::ref(handler));
    std::this_thread::sleep_for(std::chrono::seconds(1));
    
    // std::thread plotThread(rawDataPlottingLoop, std::ref(handler));
    // std::thread plotThread(processedDataPlottingLoop, std::ref(handler));

    // std::thread dataProcessorThread(dataProcessingLoop, std::ref(processor));
    std::thread dataProcessorThread(dataProcessingLoop, std::ref(handler));

    dataThread.join();
    // plotThread.join();
    dataProcessorThread.join();

    return 0;
}
