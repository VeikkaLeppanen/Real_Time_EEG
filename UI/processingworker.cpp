#include "processingworker.h"
#include <iostream>
#include <QThread>


ProcessingWorker::ProcessingWorker(dataHandler &handler, Eigen::MatrixXd &processed_data, volatile std::sig_atomic_t &processingWorkerRunning, QObject* parent)
    : QObject(parent), handler(handler), processed_data(processed_data), processingWorkerRunning(processingWorkerRunning)
{
    // Filtering
    int order = 4;
    double Fs = 5000;  // Sampling frequency
    double Fc1 = 0.33;   // Desired cutoff frequency
    double Fc2 = 135;   // Desired cutoff frequency
    int numTaps = 301;  // Length of the FIR filter
    // filterCoeffs_ = designLowPassFilter(numTaps, Fs, Fc2);
    filterCoeffs_ = designBandPassFilter(numTaps, Fs, Fc1, Fc2);
}

ProcessingWorker::~ProcessingWorker()
{ }

void ProcessingWorker::process()
{
    // Ensure this code is thread-safe and does not interfere with the GUI thread
    try {
        std::cout << "Pworker start" << '\n';
        int seq_num_tracker = 0;

        int n_eeg_channels = handler.get_channel_count();
        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_display);
        Eigen::MatrixXd EEG_filtered = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_display);

        int downsampling_factor = 10;
        int newCols = (samples_to_display + downsampling_factor - 1) / downsampling_factor;
        Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_eeg_channels, newCols);

        while(processingWorkerRunning) {
            // Initializing parameters and matrices for algorithms
            int n_eeg_channels = handler.get_channel_count();

            // Eigen::MatrixXd all_channels = handler.returnLatestDataInOrder(samples_to_display);
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_display);

            // Check how many samples skipped. REMOVE false in order to use
            // if(false && sequence_number > seq_num_tracker + 1) { 
            //     std::cout << "Number of samples skipped during processing: " << sequence_number - seq_num_tracker << '\n';
            // }

            // Filtering
            EEG_filtered = applyFIRFilterToMatrix(all_channels, filterCoeffs_);

            // Downsampling
            downsample(EEG_filtered, EEG_downsampled, downsampling_factor);

            seq_num_tracker = sequence_number;
            processed_data = EEG_downsampled;
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred: %1").arg(e.what()));
    }
}

void ProcessingWorker::process_testing()
{
    // Ensure this code is thread-safe and does not interfere with the GUI thread
    try {
        std::cout << "Pworker start" << '\n';
        int seq_num_tracker = 0;
        int count = 0;
        double total_time = 0.0;
        Eigen::MatrixXd EEG_output = Eigen::MatrixXd::Zero(127, samples_to_display);

        // Initializing parameters and matrices for algorithms
        int n_eeg_channels = handler.get_channel_count();

        // Eigen::MatrixXd all_channels = handler.returnLatestDataInOrder(samples_to_display);
        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_display);
        Eigen::MatrixXd EEG_filtered = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_display);
        
        
        while(processingWorkerRunning) {
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_display);



            // std::cout << "SeqNum: " << sequence_number << '\n';
            if (sequence_number > 1 && sequence_number % 10000 == 0 && sequence_number != seq_num_tracker) {

                std::cout << "SeqNum: " << sequence_number << '\n';
                auto start = std::chrono::high_resolution_clock::now();

                EEG_filtered = applyFIRFilterToMatrix(all_channels, filterCoeffs_);

                // Downsampling
                // int downsampling_factor = 10;
                // int newCols = (samples_to_display + downsampling_factor - 1) / downsampling_factor;
                // Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_eeg_channels, newCols);
                // downsample(EEG_filtered, EEG_downsampled, downsampling_factor);

                EEG_output.row(count) = EEG_filtered.row(0);
                std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
                std::cout << "Time taken: " << elapsed.count() << " seconds." << std::endl;
                total_time += elapsed.count();
                count++;
            }
            processed_data = EEG_filtered;

            seq_num_tracker = sequence_number;
            QThread::usleep(1); // Short sleep to prevent a tight loop, adjust as needed
        }

        if (count > 0) {
            double average_time = total_time / count;
            std::cout << "Total time taken: " << total_time << " seconds." << std::endl;
            std::cout << "Average time taken: " << average_time << " seconds." << std::endl;

            // writeMatrixToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/Betas.csv", Betas);
            writeMatrixToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/EEG_corrected.csv", EEG_output);
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred: %1").arg(e.what()));
    }
}
