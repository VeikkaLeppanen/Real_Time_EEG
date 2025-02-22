#include "preProcessingWorker.h"
#include <iostream>
#include <QThread>


void signalHandlerPrep(int signal) {
    std::cerr << "Error: Segmentation fault (signal " << signal << ")\n";
    std::exit(signal);  // Exit the program
}

preProcessingWorker::preProcessingWorker(dataHandler &handler, 
                    volatile std::sig_atomic_t &processingWorkerRunning, 
                       preprocessingParameters &prepParams_in, 
                                       QObject *parent)
    : QObject(parent), 
      handler(handler), 
      processingWorkerRunning(processingWorkerRunning)
{ 
    setPreprocessingParameters(prepParams_in);
}

preProcessingWorker::~preProcessingWorker()
{ }

void preProcessingWorker::setPreprocessingParameters(preprocessingParameters newParams) {
    std::cout << "New preprocessingParameters Parameters: " << newParams << std::endl;
    currentPrepParams = newParams;

    n_channels = handler.get_channel_count();

    samples_to_process = newParams.numberOfSamples;
    downsampling_factor = newParams.downsampling_factor;
    downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;
    delay = newParams.delay;

    print_debug("Params initialized");

    // Memory preallocation for preprocessing matrices
    all_channels = Eigen::MatrixXd::Zero(n_channels, samples_to_process);
    EEG_downsampled = Eigen::MatrixXd::Zero(n_channels, downsampled_cols);
    expCWL = Eigen::MatrixXd::Zero(n_CWL_channels_to_use * (1+2*delay), downsampled_cols);
    pinvCWL = Eigen::MatrixXd::Zero(downsampled_cols, downsampled_cols);
    EEG_corrected = Eigen::MatrixXd::Zero(n_EEG_channels_to_use, downsampled_cols);

    // Input triggers
    triggers_A = Eigen::VectorXi::Zero(samples_to_process);
    triggers_B = Eigen::VectorXi::Zero(samples_to_process);
    triggers_out = Eigen::VectorXi::Zero(samples_to_process);
    time_stamps = Eigen::VectorXd::Zero(samples_to_process);

    EEG_win_data_to_display = Eigen::MatrixXd::Zero(n_channels, samples_to_process);

    print_debug("Memory allocated");
}

void preProcessingWorker::process()
{
    std::signal(SIGSEGV, signalHandlerPrep);

    int number_of_threads = std::thread::hardware_concurrency() - 4;
    Eigen::setNbThreads(number_of_threads);
    omp_set_num_threads(number_of_threads);
    
    int eigen_threads = Eigen::nbThreads();
    std::cout << "Eigen is using " << eigen_threads << " threads." << std::endl;

    try {
        std::cout << "preProcessingWorker start" << '\n';

        print_debug("Channel names set");

        int seq_num_tracker = 0;
        int stimulation_tracker = -1;
        while(processingWorkerRunning) {
            if (processing_pause) {
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
                continue;
            }
            // Time the removeBCG function
            auto start = std::chrono::high_resolution_clock::now();
            
            print_debug("Processing start");

            int sequence_number = handler.getLatestDataAndTriggers(all_channels, triggers_A, triggers_B, triggers_out, time_stamps, samples_to_process);
            
            // Check if current sample is processed
            if (seq_num_tracker == sequence_number) {
                print_debug("Sample is already processed");
                continue;
            }

            // Set seq_num_tracker
            if (seq_num_tracker == 0) {
                seq_num_tracker = sequence_number;
                continue;
            }

            print_debug("Checks passed");

            // Downsampling
            print_debug("Downsampling");
            if (downsampling_factor > 1) downsample(all_channels, EEG_downsampled, downsampling_factor);
            else EEG_downsampled = all_channels;

            // CWL
            print_debug("Performing removeBCG");
            if (performRemoveBCG) {
            
                print_debug("Delay Embedding");
                if (delay > 0) { delayEmbed(EEG_downsampled.middleRows(n_EEG_channels_to_use, n_CWL_channels_to_use), expCWL, delay); } 
                else { expCWL = EEG_downsampled.middleRows(n_EEG_channels_to_use, n_CWL_channels_to_use); }

                print_debug("removeBCG");
                removeBCG(EEG_downsampled.topRows(n_EEG_channels_to_use), expCWL, pinvCWL, EEG_corrected);

                // Update the EEG window graph data
                if (EEG_downsampled.cols() != EEG_win_data_to_display.cols()) EEG_win_data_to_display.resize(EEG_win_data_to_display.rows(), EEG_downsampled.cols());

                EEG_win_data_to_display.topRows(EEG_corrected.rows()) = EEG_corrected;
                EEG_win_data_to_display.bottomRows(EEG_downsampled.rows() - n_EEG_channels_to_use) = EEG_downsampled.bottomRows(EEG_downsampled.rows() - n_EEG_channels_to_use);


            } else {
                EEG_corrected = EEG_downsampled.topRows(n_EEG_channels_to_use);
                EEG_win_data_to_display = EEG_downsampled;
            }

            int cols_to_save = sequence_number - seq_num_tracker;
            if (cols_to_save > 0 && cols_to_save <= EEG_corrected.cols()) {
                // Eigen::MatrixXd output_matrix = EEG_corrected.rightCols(cols_to_save);
                emit savePreprocessingOutput(EEG_corrected.rightCols(cols_to_save));
            } else {
                print_debug("Skipping save - invalid column parameters");
            }

            seq_num_tracker = sequence_number;

            emit preprocessingOutputReady(EEG_corrected, triggers_A, triggers_B, triggers_out, time_stamps, samples_to_process, sequence_number);

            emit updateEEGDisplayedData(EEG_win_data_to_display, triggers_A, triggers_B, time_stamps, handler.getChannelNames());

            // Save C3
            // if (save_seqnum_tracker + 10000 < sequence_number) {
            //     C3_save.row(save_index) = EEG_corrected.row(0);
            //     save_index++;
            //     save_seqnum_tracker = sequence_number;
            // }
            
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;  // Explicitly specify duration type
            if(sequence_number > 100000) {
                total_bcg_time += duration;
                min_bcg_time = std::min(min_bcg_time, duration);
                max_bcg_time = std::max(max_bcg_time, duration);
                bcg_call_count++;
            }
        }

        // Print timing statistics when exiting the loop
        if (bcg_call_count > 0) {
            double avg_time = total_bcg_time.count() / bcg_call_count;
            double min_time = min_bcg_time.count();
            double max_time = max_bcg_time.count();
            
            std::cout << "removeBCG statistics:\n"
                      << "Total calls: " << bcg_call_count << "\n"
                      << "Average time: " << avg_time * 1000 << " ms\n"
                      << "Minimum time: " << min_time * 1000 << " ms\n"
                      << "Maximum time: " << max_time * 1000 << " ms\n";
        }

        // writeMatrixdToCSV("C3_save.csv", C3_save);

        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in processingworker process function: %1").arg(e.what()));
    }
}





