#include <vector>
#include <thread>
#include <mutex>
#include <omp.h>
#include <string>
#include <csignal>
#include <chrono>

#include "dataProcessor/dataProcessor.h"
#include "dataHandler/dataHandler.h"
#include "eeg_bridge/eeg_bridge.h"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

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

// Loop for matplotlibcpp graph visualization
void rawDataPlottingLoop(dataHandler &handler) {
    
    // Wait for the handler initialization with measurementStartPackage
    while (!handler.isReady()) {
        // Optionally, sleep for a short duration to avoid busy waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    int datapoints = handler.get_buffer_capacity();
    int dataPoints_downsampled = (datapoints + DOWNSAMPLING_FACTOR - 1) / DOWNSAMPLING_FACTOR;

    std::vector<int> channels_to_display = {0, 1, 2};

    plt::ion(); // Enable interactive mode
    plt::figure_size(1920, 100 * (channels_to_display.size() + 2));

    // X values of the plot
    Eigen::VectorXd xVec = Eigen::VectorXd::LinSpaced(dataPoints_downsampled, handler.get_buffer_length_in_seconds(), 0);
    std::vector<double> x(xVec.data(), xVec.data() + xVec.size());

    // Y values of the plot
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(dataPoints_downsampled);
    Eigen::MatrixXd downSampledData = Eigen::MatrixXd::Zero(handler.get_channel_count(), dataPoints_downsampled);

    while (!signal_received) {

        plt::clf();

        int graph_ind = 0;
        for (int channel_index : channels_to_display) {

            yVec = handler.getChannelDataInOrder(channel_index, DOWNSAMPLING_FACTOR);

            // std::cout << "yVec: "; 
            // for (int abc = 0; abc < yVec.size(); abc++) {
            //     std::cout << yVec[abc] << ' ';
            // }
            // std::cout << '\n';

            // std::cout << yVec.size() << ' ' << x.size() << '\n';
            plt::subplot(channels_to_display.size() + 2, 1, graph_ind + 1);
            plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
            plt::title("Channel " + std::to_string(channel_index + 1)); // Optional: Add title to each subplot
            graph_ind++;
        }
        
        // Timestamp plotting
        // yVec = handler.getTimeStampsInOrder(DOWNSAMPLING_FACTOR);
        // plt::subplot(channels_to_display.size() + 2, 1, graph_ind + 1);
        // plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
        // plt::title("Time stamps");
        // graph_ind++;

        // Trigger plotting
        yVec = handler.getTriggersInOrder(DOWNSAMPLING_FACTOR);
        plt::subplot(channels_to_display.size() + 1, 1, graph_ind + 1);
        plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
        plt::title("Triggers");
        
        plt::pause(0.01); // Pause for a short period to allow the plot to update
    }

    plt::close();
    std::cout << "Exiting plottingLoop" << '\n';
}

















// DATA PROCESSING AND PROCESSED DATA PLOTTING

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

Eigen::MatrixXd delayEmbed(const Eigen::MatrixXd& CWL, int delay) {
    int rows = CWL.rows();
    int cols = CWL.cols();
    Eigen::MatrixXd expandedCWL(rows, cols + 2 * delay);

    // Pad columns with zeros based on delay
    expandedCWL.block(0, delay, rows, cols) = CWL;

    // Fill in zeros for padding
    expandedCWL.block(0, 0, rows, delay).setZero();
    expandedCWL.block(0, cols + delay, rows, delay).setZero();

    return expandedCWL;
}



void removeBCG(const Eigen::MatrixXd& EEG, Eigen::MatrixXd& CWL, Eigen::MatrixXd& EEG_corrected, int delay) {

    if (delay > 0) {
        // Flip and perform delay embedding if delay is positive
        CWL = delayEmbed(CWL.colwise().reverse(), 1 + 2 * delay).colwise().reverse();
    }
    
    
    // EEG_corrected = EEG - (EEG * (CWL.completeOrthogonalDecomposition().pseudoInverse() * CWL));

    // Eigen::MatrixXd squareCWL = CWL.completeOrthogonalDecomposition().pseudoInverse() * CWL;


    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(CWL, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // Eigen::VectorXd singular_values = svd.singularValues();
    // Eigen::MatrixXd singular_values_inv = Eigen::MatrixXd::Zero(CWL.cols(), CWL.rows());

    // // Invert the singular values
    // for (int i = 0; i < singular_values.size(); ++i) {
    //     if (singular_values(i) > 1e-10) {  // Threshold to handle very small values
    //         singular_values_inv(i, i) = 1.0 / singular_values(i);
    //     }
    // }

    // Eigen::MatrixXd A_pinv = svd.matrixV() * singular_values_inv * svd.matrixU().transpose();


    // Calculate the pseudoinverse using the formula A^+ = (A^T A)^{-1} A^T
    Eigen::MatrixXd CWLtCWL = CWL.transpose() * CWL; // Compute A^T A
    Eigen::MatrixXd CWLtCWL_inv = CWLtCWL.inverse(); // Compute (A^T A)^{-1}
    Eigen::MatrixXd CWL_pinv = CWLtCWL_inv * CWL.transpose(); // Compute pseudoinverse


    Eigen::MatrixXd squareCWL = CWL_pinv * CWL;

    Eigen::MatrixXd EEG_fits = Eigen::MatrixXd::Zero(EEG.rows(), EEG.cols());
    for(int i = 0; i < EEG.rows(); i++) {
        EEG_fits.row(i) = EEG.row(i) * squareCWL;
    }

    EEG_corrected = EEG - EEG_fits;
    // std::cout << EEG_corrected.row(0) << '\n';
}

// Loop for the data processing
void dataProcessingLoop(dataHandler &handler) {
    
    int n_eeg_channels = 5;
    int n_cwl_channels = 7;
    int n_samples = 5000;
    Eigen::MatrixXd EEG = Eigen::MatrixXd::Zero(n_eeg_channels, n_samples);
    Eigen::MatrixXd CWL = Eigen::MatrixXd::Zero(n_cwl_channels, n_samples);

    Eigen::MatrixXd EEG_filtered = Eigen::MatrixXd::Zero(n_eeg_channels, n_samples);
    Eigen::MatrixXd CWL_filtered = Eigen::MatrixXd::Zero(n_cwl_channels, n_samples);

    int downsampling_factor = 10;
    int newCols = (n_samples + downsampling_factor - 1) / downsampling_factor;
    Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_eeg_channels, newCols);
    Eigen::MatrixXd CWL_downsampled = Eigen::MatrixXd::Zero(n_cwl_channels, newCols);

    int delay = 0;
    Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_eeg_channels, newCols);

    while (!handler.isReady()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    while (!signal_received) {
        auto start = std::chrono::high_resolution_clock::now();

        // Copy data from wanted EEG channels and CWL channels
        EEG = handler.getBlockChannelDataInOrder(0, n_eeg_channels, n_samples);
        CWL = handler.getBlockChannelDataInOrder(5, n_cwl_channels, n_samples);


        // LOWPASS filter 120Hz   BANDPASS 0.33-125Hz FIR
        EEG_filtered = EEG;
        CWL_filtered = CWL;

        // DONWSAMPLING
        downsample(EEG_filtered, EEG_downsampled, downsampling_factor);
        downsample(CWL_filtered, CWL_downsampled, downsampling_factor);

        // auto start = std::chrono::high_resolution_clock::now();
        // // removeBCG
        // removeBCG(EEG_downsampled, CWL_downsampled, EEG_corrected, delay);

        std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
        std::cout << "Time taken: " << elapsed.count() << " seconds." << std::endl;
    }

    std::cout << "Exiting dataProcessingLoop" << '\n';
}

// Loop for matplotlibcpp graph visualization
void processedDataPlottingLoop(dataHandler &handler) {
    
    // Wait for the handler initialization with measurementStartPackage
    while (!handler.isReady()) {
        // Optionally, sleep for a short duration to avoid busy waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    int datapoints = handler.get_buffer_capacity();
    int dataPoints_downsampled = (datapoints + DOWNSAMPLING_FACTOR - 1) / DOWNSAMPLING_FACTOR;

    std::vector<int> channels_to_display = {0, 1, 2};

    plt::ion(); // Enable interactive mode
    plt::figure_size(1920, 100 * (channels_to_display.size() + 2));

    // X values of the plot
    Eigen::VectorXd xVec = Eigen::VectorXd::LinSpaced(dataPoints_downsampled, handler.get_buffer_length_in_seconds(), 0);
    std::vector<double> x(xVec.data(), xVec.data() + xVec.size());

    // Y values of the plot
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(dataPoints_downsampled);
    Eigen::MatrixXd downSampledData = Eigen::MatrixXd::Zero(handler.get_channel_count(), dataPoints_downsampled);

    while (!signal_received) {

        plt::clf();

        int graph_ind = 0;
        for (int channel_index : channels_to_display) {

            yVec = handler.getChannelDataInOrder(channel_index, DOWNSAMPLING_FACTOR);

            plt::subplot(channels_to_display.size() + 2, 1, graph_ind + 1);
            plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
            plt::title("Channel " + std::to_string(channel_index + 1));
            graph_ind++;
        }
        
        plt::pause(0.01); // Pause for a short period to allow the plot to update
    }

    plt::close();
    std::cout << "Exiting plottingLoop" << '\n';
}






















int main() {
    std::signal(SIGINT, signal_handler);
    
    // size_t thread_count = 6;
    // dataProcessor processor(thread_count);
    // dataHandler handler(processor);
    dataHandler handler;

    // std::thread dataThread(dataSimulationLoop, std::ref(handler));
    std::thread dataThread(dataAcquisitionLoop, std::ref(handler));
    std::this_thread::sleep_for(std::chrono::seconds(1));
    
    std::thread plotThread(rawDataPlottingLoop, std::ref(handler));
    // std::thread plotThread(processedDataPlottingLoop, std::ref(handler));

    // std::thread dataProcessorThread(dataProcessingLoop, std::ref(processor));
    std::thread dataProcessorThread(dataProcessingLoop, std::ref(handler));

    dataThread.join();
    plotThread.join();
    dataProcessorThread.join();

    return 0;
}

















