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
    omp_set_num_threads(1);
    
    int eigen_threads = Eigen::nbThreads();
    std::cout << "Eigen is using " << eigen_threads << " threads." << std::endl;

    try {
        std::cout << "preProcessingWorker start" << '\n';

        print_debug("Channel names set");

        // int index = 0;
        // std::vector<int> seqNum_list;
        // Eigen::MatrixXd EEG_save_1 = Eigen::MatrixXd::Zero(200, downsampled_cols);
        // Eigen::MatrixXd EEG_save_2 = Eigen::MatrixXd::Zero(200, downsampled_cols);

        int seq_num_tracker = 0;
        int stimulation_tracker = -1;
        while(processingWorkerRunning) {
            if (processing_pause) {
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
                continue;
            }
            
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
            print_debug("Downsampled");
            if (downsampling_factor > 1) downsample(all_channels, EEG_downsampled, downsampling_factor);
            else EEG_downsampled = all_channels;

            // CWL
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
                
                // if (sequence_number % 10000 == 0) {
                //     EEG_save_1.row(index) = EEG_corrected.row(0);
                //     EEG_save_2.row(index) = EEG_corrected.row(0) - EEG_corrected.bottomRows(4).colwise().mean();
                //     seqNum_list.push_back(sequence_number);
                //     std::cout << "index saved: " << index << " SeqNum: " << sequence_number << '\n';
                //     index++;
                // }

                // std::this_thread::sleep_for(std::chrono::milliseconds(3));
            }

            seq_num_tracker = sequence_number;

            // handler.setPreprocessingOutput(EEG_corrected, triggers_A, triggers_B, triggers_out, time_stamps, samples_to_process, sequence_number);
            emit preprocessingOutputReady(EEG_corrected, triggers_A, triggers_B, triggers_out, time_stamps, samples_to_process, sequence_number);

            emit updateEEGDisplayedData(EEG_win_data_to_display, triggers_A, triggers_B, time_stamps, handler.getChannelNames());
        }

        // writeMatrixdToCSV("data_reference.csv", EEG_save_1);
        // writeMatrixdToCSV("data_reference_spatial.csv", EEG_save_2);
        // writeMatrixiToCSV("seqNum_list.csv", vectorToColumnMatrixi(seqNum_list));

        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in processingworker process function: %1").arg(e.what()));
    }
}





