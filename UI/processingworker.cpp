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
    // Set the current thread to use SCHED_RR
    pthread_t this_thread = pthread_self();
    struct sched_param params;
    params.sched_priority = sched_get_priority_max(SCHED_RR);
    
    if (pthread_setschedparam(this_thread, SCHED_RR, &params) != 0) {
        qDebug("Failed to set thread to real-time");
    } else {
        qDebug("Thread set to real-time successfully");
    }
    
    try {
        std::cout << "Pworker start" << '\n';
        int seq_num_tracker = 0;

        //  downsampling
        int downsampling_factor = 10;
        int downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // removeBCG
        int delay = 6;

        // phase estimate
        Eigen::VectorXd LSFIR_coeffs_1;
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_0_80Hz(LSFIR_coeffs_1);
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);

        size_t edge = 35;
        size_t modelOrder = 15;
        size_t hilbertWinLength = 64;
        size_t estimationLength = edge + std::ceil(hilbertWinLength / 2);

        // Initializing parameters and matrices for algorithms
        int n_eeg_channels = handler.get_channel_count();

        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_process);
        Eigen::MatrixXd EEG_filter1 = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_process);
        Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_eeg_channels, downsampled_cols);
        Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_eeg_channels, downsampled_cols);
        Eigen::VectorXd EEG_spatial = Eigen::VectorXd::Zero(downsampled_cols);
        Eigen::VectorXd EEG_filter2 = Eigen::VectorXd::Zero(downsampled_cols);
        std::vector<double> EEG_predicted(estimationLength, 0.0);
        std::vector<std::complex<double>> EEG_hilbert(estimationLength, std::complex<double>(0.0, 0.0));
        Eigen::VectorXd phaseAngles(estimationLength);

        Eigen::MatrixXd Data_to_display(1, estimationLength);

        
        while(processingWorkerRunning) {
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_process);

            std::cout << "Samples skipped: " << sequence_number - seq_num_tracker << '\n';

            // Filtering
            // EEG_filter1 = applyFIRFilterToMatrix(all_channels * 10, filterCoeffs_);                        // Testing scaling here
            EEG_filter1 = zeroPhaseLSFIRMatrix(all_channels * 10, LSFIR_coeffs_1);

            // Downsampling
            downsample(EEG_filter1, EEG_downsampled, downsampling_factor);

            // CWL
            removeBCG(EEG_downsampled.middleRows(0, 5), EEG_downsampled.middleRows(5, 7), EEG_corrected, delay);
            EEG_spatial = EEG_corrected.row(0) - EEG_corrected.bottomRows(4).colwise().mean();

            // Phase estimate
            EEG_filter2 = zeroPhaseLSFIR(EEG_spatial, LSFIR_coeffs_2);

            // Hilbert transform
            EEG_predicted = fitAndPredictAR_LeastSquares(EEG_filter2, modelOrder, estimationLength);
            EEG_hilbert = hilbertTransform(EEG_predicted);
            for (std::size_t i = 0; i < estimationLength; ++i) {
                phaseAngles[i] = std::arg(EEG_hilbert[i]);
            }

            // Save displayed data
            Data_to_display.row(0) = phaseAngles;
            processed_data = Data_to_display;

            seq_num_tracker = sequence_number;
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in process_test function: %1").arg(e.what()));
    }
}

void ProcessingWorker::process_testing()
{
    // Set the current thread to use SCHED_RR
    pthread_t this_thread = pthread_self();
    struct sched_param params;
    params.sched_priority = sched_get_priority_max(SCHED_RR);
    
    if (pthread_setschedparam(this_thread, SCHED_RR, &params) != 0) {
        qDebug("Failed to set thread to real-time");
    } else {
        qDebug("Thread set to real-time successfully");
    }

    // Ensure this code is thread-safe and does not interfere with the GUI thread
    try {
        std::cout << "Pworker start" << '\n';
        int seq_num_tracker = 0;
        int count = 0;
        double filtering_total_time = 0.0;
        double downsampling_total_time = 0.0;
        double removeBCG_total_time = 0.0;
        double pEstFilt_total_time = 0.0;
        double phaseEstimate_total_time = 0.0;
        double hilbert_total_time = 0.0;
        double total_time = 0.0;

        //  downsampling
        int downsampling_factor = 10;
        int downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // removeBCG
        int delay = 6;

        // filters
        Eigen::VectorXd LSFIR_coeffs_1;
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_0_80Hz(LSFIR_coeffs_1);
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);

        Eigen::VectorXd a, b;
        getBWCoeffs_8Hz_12Hz(a, b);

        size_t edge = 35;
        int edge_cut_cols = downsampled_cols - 2 * edge;
        size_t modelOrder = 15;
        size_t hilbertWinLength = 64;
        size_t estimationLength = edge + std::ceil(hilbertWinLength / 2);

        // Initializing parameters and matrices for algorithms
        int n_eeg_channels = handler.get_channel_count();

        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_process);
        Eigen::MatrixXd EEG_filter1 = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_process);
        Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_eeg_channels, downsampled_cols);
        Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_eeg_channels, downsampled_cols);
        Eigen::VectorXd EEG_spatial = Eigen::VectorXd::Zero(downsampled_cols);
        Eigen::VectorXd EEG_filter2 = Eigen::VectorXd::Zero(downsampled_cols);
        std::vector<double> EEG_predicted(estimationLength, 0.0);
        std::vector<std::complex<double>> EEG_hilbert(estimationLength, std::complex<double>(0.0, 0.0));
        Eigen::VectorXd phaseAngles(estimationLength);

        Eigen::MatrixXd CSV_save = Eigen::MatrixXd::Zero(164, downsampled_cols);


        Eigen::MatrixXd Data_to_display(3, downsampled_cols);

        
        while(processingWorkerRunning) {
            
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_process);

            // std::cout << "SeqNum: " << sequence_number << '\n';
            if (sequence_number > 1 && sequence_number % 10000 == 0 && sequence_number != seq_num_tracker) {

                std::cout << "SeqNum: " << sequence_number << " count: " << count << '\n';
                auto start = std::chrono::high_resolution_clock::now();

                // Filtering
                // EEG_filter1 = applyFIRFilterToMatrix(all_channels * 10, filterCoeffs_);                        // Testing scaling here
                EEG_filter1 = zeroPhaseLSFIRMatrix(all_channels * 10, LSFIR_coeffs_1);

                std::cout << "After first filter!" << '\n';                
                std::chrono::duration<double> filtering_time = std::chrono::high_resolution_clock::now() - start;
                filtering_total_time += filtering_time.count();

                // Downsampling
                downsample(EEG_filter1, EEG_downsampled, downsampling_factor);

                std::chrono::duration<double> downsampling_time = std::chrono::high_resolution_clock::now() - filtering_time - start;
                downsampling_total_time += downsampling_time.count();

                // CWL
                removeBCG(EEG_downsampled.middleRows(0, 5), EEG_downsampled.middleRows(5, 7), EEG_corrected, delay);
                EEG_spatial = EEG_corrected.row(0) - EEG_corrected.bottomRows(4).colwise().mean();

                std::chrono::duration<double> removeBCG_time = std::chrono::high_resolution_clock::now() - downsampling_time - filtering_time - start;
                removeBCG_total_time += removeBCG_time.count();

                // Second filtering
                EEG_filter2 = zeroPhaseLSFIR(EEG_spatial, LSFIR_coeffs_2);
                // EEG_filter2 = zeroPhaseBW(EEG_spatial, a, b);

                std::chrono::duration<double> pEstFilt_time = std::chrono::high_resolution_clock::now() - removeBCG_time - downsampling_time - filtering_time - start;
                pEstFilt_total_time += pEstFilt_time.count();

                // Phase estimate
                EEG_predicted = fitAndPredictAR_LeastSquares(EEG_filter2.segment(edge, edge_cut_cols), modelOrder, estimationLength);

                std::chrono::duration<double> phaseEstimate_time = std::chrono::high_resolution_clock::now() - pEstFilt_time - removeBCG_time - downsampling_time - filtering_time - start;
                phaseEstimate_total_time += phaseEstimate_time.count();

                // Hilbert transform
                EEG_hilbert = hilbertTransform(EEG_predicted);
                for (std::size_t i = 0; i < estimationLength; ++i) {
                    phaseAngles[i] = std::arg(EEG_hilbert[i]);
                }

                std::chrono::duration<double> hilbert_time = std::chrono::high_resolution_clock::now() - phaseEstimate_time - pEstFilt_time - removeBCG_time - downsampling_time - filtering_time - start;
                hilbert_total_time += hilbert_time.count();

                // Save csv and displayed data
                CSV_save.row(count) = EEG_filter2;
                Data_to_display.row(0) = EEG_downsampled.row(0);
                Data_to_display.row(1) = EEG_spatial;
                Data_to_display.row(2) = EEG_filter2;
                processed_data = Data_to_display;
                
                std::chrono::duration<double> total_elapsed = std::chrono::high_resolution_clock::now() - start;
                std::cout << "Time taken: " << total_elapsed.count() << " seconds." << std::endl;
                total_time += total_elapsed.count();
                count++;
            }

            seq_num_tracker = sequence_number;
            QThread::usleep(1); // Short sleep to prevent a tight loop, adjust as needed
        }

        if (count > 0) {
            double filtering_average_time = filtering_total_time / count;
            double downsampling_average_time = downsampling_total_time / count;
            double removeBCG_average_time = removeBCG_total_time / count;
            double pEstFilt_average_time = pEstFilt_total_time / count;
            double phaseEstimate_average_time = phaseEstimate_total_time / count;
            double hilbert_average_time = hilbert_total_time / count;
            double average_time = total_time / count;
            std::cout << "Total time taken: " << total_time << " seconds." << std::endl;
            std::cout << "Average filtering time taken: " << filtering_average_time << " seconds." << std::endl;
            std::cout << "Average downsampling time taken: " << downsampling_average_time << " seconds." << std::endl;
            std::cout << "Average removeBCG time taken: " << removeBCG_average_time << " seconds." << std::endl;
            std::cout << "Average second filtering time taken: " << pEstFilt_average_time << " seconds." << std::endl;
            std::cout << "Average phase estimate time taken: " << phaseEstimate_average_time << " seconds." << std::endl;
            std::cout << "Average hilbert time taken: " << hilbert_average_time << " seconds." << std::endl;
            std::cout << "Average total time taken: " << average_time << " seconds." << std::endl;

            // writeMatrixToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/Betas.csv", Betas);
            writeMatrixToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/EEG_corrected.csv", CSV_save);
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in process_test function: %1").arg(e.what()));
    }
}
