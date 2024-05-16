#include "processingworker.h"
#include <iostream>

ProcessingWorker::ProcessingWorker(dataHandler &handler, Eigen::MatrixXd &processed_data, volatile std::sig_atomic_t &processingWorkerRunning, QObject* parent)
    : QObject(parent), handler(handler), processed_data(processed_data), processingWorkerRunning(processingWorkerRunning)
{
    // Filtering
    int order = 4;
    double Fs = 5000;  // Sampling frequency
    double Fc1 = 8;   // Desired cutoff frequency
    double Fc2 = 120;   // Desired cutoff frequency
    int numTaps = 51;  // Length of the FIR filter
    // filterCoeffs_ = designLowPassFilter(numTaps, Fs, Fc2);
    filterCoeffs_ = designBandPassFilter(numTaps, Fs, Fc1, Fc2);
    // computeButterworthCoefficients(order, Fs, Fc, b, a);
}

ProcessingWorker::~ProcessingWorker()
{ }

void ProcessingWorker::process()
{
    // Ensure this code is thread-safe and does not interfere with the GUI thread
    try {
        std::cout << "Pworker start" << '\n';

        while(processingWorkerRunning) {
        // Initializing parameters and matrices for algorithms
        int n_eeg_channels = handler.get_channel_count();

        // Eigen::MatrixXd all_channels = handler.returnLatestDataInOrder(samples_to_display);
        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_display);
        int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_display);

        // std::cout << "SeqNum: " << sequence_number << '\n';
        // if (sequence_number > 1 && sequence_number % 10000 == 0) {
        //     std::cout << "SeqNum: " << sequence_number << '\n';

            // Filtering
            Eigen::MatrixXd EEG_filtered = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_display);
            auto start = std::chrono::high_resolution_clock::now();

            // EEG_filtered = applyFIRFilterToMatrix(all_channels, filterCoeffs_);
            EEG_filtered = applyBandFIRFilterToMatrix(all_channels, filterCoeffs_);
            // EEG_filtered = applyIIRFilterToMatrix(all_channels, b, a);

            // Downsampling
            // int downsampling_factor = 10;
            // int newCols = (samples_to_display + downsampling_factor - 1) / downsampling_factor;
            // Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_eeg_channels, newCols);
            // downsample(EEG_filtered, EEG_downsampled, downsampling_factor);

            // std::cout << EEG_output.cols() << ' ' << EEG_filtered.cols() << '\n';
            // EEG_output.row(count) = EEG_filtered.row(0);
            // std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
            // std::cout << "Time taken: " << elapsed.count() << " seconds." << std::endl;
            // total_time += elapsed.count();
            // count++;
            processed_data = EEG_filtered;
        // }
        }
            
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred: %1").arg(e.what()));
    }
}
