#include "processingworker.h"
#include <iostream>
#include <QThread>


void signalHandler(int signal) {
    std::cerr << "Error: Segmentation fault (signal " << signal << ")\n";
    std::exit(signal);  // Exit the program
}

ProcessingWorker::ProcessingWorker(dataHandler &handler, 
                               Eigen::MatrixXd &processed_data, 
                    volatile std::sig_atomic_t &processingWorkerRunning, 
                       preprocessingParameters &prepParams, 
                       phaseEstimateParameters &phaseEstParams, 
                                       QObject *parent)
    : QObject(parent), 
      handler(handler), 
      processed_data(processed_data), 
      processingWorkerRunning(processingWorkerRunning), 
      prepParams(prepParams),
      phaseEstParams(phaseEstParams)
{ }

ProcessingWorker::~ProcessingWorker()
{ }

void ProcessingWorker::process()
{
    std::signal(SIGSEGV, signalHandler);
    omp_set_num_threads(10);

    // Eigen::setNbThreads(std::thread::hardware_concurrency());

    try {
        std::cout << "ProcessingWorker start" << '\n';
        std::cout << "prepParams: " << prepParams << '\n';
        std::cout << "phaseEstParams: " << phaseEstParams << '\n';

        // Number of samples to use for processing
        int samples_to_process = prepParams.numberOfSamples;

        //  downsampling
        int downsampling_factor = prepParams.downsampling_factor;
        downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // removeBCG
        int delay = prepParams.delay;
        int n_EEG_channels_to_use = 5;      
        int n_CWL_channels_to_use = 7;

        // FIR filters
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);

        // phase estimate
        setPhaseEstimateParameters(phaseEstParams);

        int n_channels = handler.get_channel_count();
        
        print_debug("Params initialized");

        // Memory preallocation for preprocessing matrices
        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_channels, samples_to_process);
        Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_channels, downsampled_cols);
        Eigen::MatrixXd expCWL = Eigen::MatrixXd::Zero(n_CWL_channels_to_use * (1+2*delay), downsampled_cols);
        Eigen::MatrixXd pinvCWL = Eigen::MatrixXd::Zero(downsampled_cols, downsampled_cols);
        Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_EEG_channels_to_use, downsampled_cols);
        Eigen::VectorXd EEG_spatial = Eigen::VectorXd::Zero(downsampled_cols);

        // Input triggers
        Eigen::VectorXi triggers_A = Eigen::VectorXi::Zero(samples_to_process);
        Eigen::VectorXi triggers_B = Eigen::VectorXi::Zero(samples_to_process);
        Eigen::VectorXi triggers_out = Eigen::VectorXi::Zero(samples_to_process);

        print_debug("Memory allocated");

        // Set names for each channel in Data_to_display
        std::vector<std::string> EEG_channel_names;
        std::vector<std::string> PhaseEst_channel_names;
        print_debug("Channel names set");

        int seq_num_tracker = 0;
        int stimulation_tracker = -1;
        while(processingWorkerRunning) {
            print_debug("Processing start");
            if (!performPreprocessing) continue;
    
            PhaseEst_channel_names.clear();
            EEG_channel_names = handler.getChannelNames();

            if (EEG_channel_names.size() >= n_EEG_channels_to_use) { 
                std::vector<std::string> EEG_spatial_channel_names(EEG_channel_names.begin(), EEG_channel_names.begin() + n_EEG_channels_to_use);
                emit updateSpatialChannelNames(EEG_spatial_channel_names);
            }

            int sequence_number = handler.getLatestDataAndTriggers(all_channels, triggers_A, triggers_B, triggers_out, samples_to_process);
            
            // Check if current sample is processed
            if (seq_num_tracker == sequence_number) continue;

            // Set seq_num_tracker
            if (seq_num_tracker == 0) {
                seq_num_tracker = sequence_number;
                continue;
            }

            print_debug("Checks passed");

            // Downsampling
            print_debug("Downsampled");
            downsample(all_channels, EEG_downsampled, downsampling_factor);

            // CWL
            if (performRemoveBCG) {
                print_debug("Delay Embedding");
                if (delay > 0) { delayEmbed(EEG_downsampled.middleRows(n_EEG_channels_to_use, n_CWL_channels_to_use), expCWL, delay); } 
                else { expCWL = EEG_downsampled.middleRows(n_EEG_channels_to_use, n_CWL_channels_to_use); }

                print_debug("removeBCG");
                removeBCG(EEG_downsampled.topRows(n_EEG_channels_to_use), expCWL, pinvCWL, EEG_corrected);
            } 
            else {
                EEG_corrected = EEG_downsampled.topRows(n_EEG_channels_to_use);
            }

            // Update the EEG window graph data
            switch (display_state) 
            {
                case displayState::RAW:
                    emit updateEEGDisplayedData(all_channels, triggers_A, triggers_B);
                    emit updateEEGwindowNames(EEG_channel_names);
                    break;
                case displayState::DOWNSAMPLED:
                    emit updateEEGDisplayedData(EEG_downsampled, triggers_A, triggers_B);
                    emit updateEEGwindowNames(EEG_channel_names);
                    break;
                case displayState::CWL:
                    emit updateEEGDisplayedData(EEG_corrected, triggers_A, triggers_B);
                    emit updateEEGwindowNames(EEG_channel_names);
                    break;
            }

            if (!performPhaseEstimation) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                continue;
            }

            emit sendNumSamples(samples_to_process);

            print_debug("Spatial filtering");
            if (spatial_channel_index >= 0 && spatial_channel_index < EEG_corrected.rows()) {
                
                Eigen::VectorXd sum_of_rows = Eigen::VectorXd::Zero(downsampled_cols);
                int outer_channel_index = 0;
                for (int i = 0; i < EEG_corrected.rows(); i++) {
                    if (i == spatial_channel_index) continue;
                        
                    if (outerElectrodeCheckStates_[outer_channel_index]) {
                            sum_of_rows += EEG_corrected.row(i);
                    }
                    outer_channel_index++;
                }
                
                int trueCount = std::count(outerElectrodeCheckStates_.begin(), outerElectrodeCheckStates_.end(), true);
                if (trueCount > 0) {
                    EEG_spatial = EEG_corrected.row(spatial_channel_index) - (sum_of_rows / trueCount).transpose();
                } else {
                    EEG_spatial = EEG_corrected.row(spatial_channel_index);
                }

            } else {
                // Handle the case where the index is out of bounds
                std::cerr << "Error: Row index is out of bounds." << std::endl;
            }
            
            // SNR check
            if (performSNRcheck) {
                print_debug("SNR check");
                double SNR = calculateSNR(EEG_spatial, 64, 500, 500.0, 10.0, 2.0);
                if (SNR < 1.0) {
                    std::cout << "SNR too small: " << SNR << '\n';
                    seq_num_tracker = sequence_number;
                    continue;
                }
            }

            // Filter2
            print_debug("Second filtering");
            if (performFiltering) {
                EEG_filter2 = zeroPhaseLSFIR(EEG_spatial.tail(filter2_length), LSFIR_coeffs_2);
            } else {
                EEG_filter2 = EEG_spatial.tail(filter2_length);
            }

            // Phase estimation
            print_debug("AR predicted");
            if (performEstimation) {
                EEG_predicted = fitAndPredictAR_YuleWalker_V2(EEG_filter2.segment(phaseEstParams.edge, edge_cut_cols), phaseEstParams.modelOrder, estimationLength);
            }
            
            // Hilbert transform
            print_debug("Hilbert transform");
            if (performHilbertTransform) {
                EEG_hilbert = hilbertTransform(EEG_predicted);
            }

            // Trigger phase targeting
            print_debug("Phase targeting");
            if (performPhaseTargeting) {
                int trigger_seqNum = findTargetPhase(EEG_hilbert, phaseAngles, sequence_number, downsampling_factor, phaseEstParams.edge, phaseEstParams.phase_shift, phaseEstParams.stimulation_target);
                if (trigger_seqNum) { 
                    handler.insertTrigger(trigger_seqNum); 
                    if (stimulation_tracker < 0) stimulation_tracker++;
                }
            }

            Data_to_display.topLeftCorner(5, downsampled_cols) = EEG_corrected;
            Data_to_display.row(5).head(downsampled_cols) = EEG_spatial;
            Data_to_display.row(6).segment(downsampled_cols - filter2_length, filter2_length - phaseEstParams.edge) = EEG_filter2.head(filter2_length - phaseEstParams.edge);

            Data_to_display.row(6).tail(estimationLength) = Eigen::Map<Eigen::VectorXd>(EEG_predicted.data(), EEG_predicted.size());
            
            int phase_length = 32;
            int phase_start = phaseEstParams.edge - phase_length / 2;
            Data_to_display.row(7).segment(downsampled_cols - phase_length / 2, phase_length) = phaseAngles.segment(phase_start, phase_length);

            if (performPhaseDifference) {
                int numSkippedSamples = (sequence_number - last_phase_seqnum) / downsampling_factor;              // FIGURE OUT A BETTER WAY TO HANDLE DOWNSAMPLING
                if (last_phase_seqnum == -1) {
                    last_phase = phaseAngles(phaseEstParams.edge);
                    last_phase_seqnum = sequence_number;
                } else if (numSkippedSamples > 35) {
                    double difference;
                    switch (phaseErrorType) {
                        case 0:
                            phase_diff_hilbert = hilbertTransform(EEG_filter2);
                            difference = ang_diff(last_phase, std::arg(phase_diff_hilbert[filter2_length - numSkippedSamples]));

                            // Update the phase difference
                            phaseDifference(phaseDifference_current_index) = difference;

                            phaseDifference_current_index = (phaseDifference_current_index + 1) % phaseDifference.size();
                            break;
                        case 1:
                            phase_diff_hilbert = hilbertTransform(EEG_filter2);
                            difference = ang_diff(last_phase, std::arg(phase_diff_hilbert[filter2_length - (numSkippedSamples)]));

                            // Update the phase difference
                            if (phaseDifference_current_index + numSkippedSamples < phaseDifference.size()) {
                                // Safe to modify the segment without wrapping around
                                phaseDifference.segment(phaseDifference_current_index, numSkippedSamples).setConstant(difference);
                            } else {
                                // Calculate the number of elements that fit until the end of the vector
                                int fitToEnd = phaseDifference.size() - phaseDifference_current_index;

                                // Update the part that fits
                                phaseDifference.segment(phaseDifference_current_index, fitToEnd).setConstant(difference);

                                // Calculate the remaining elements that need to wrap around to the beginning
                                int overflow = numSkippedSamples - fitToEnd;

                                // Update the wrapped around part
                                if (overflow > 0) {
                                    phaseDifference.head(overflow).setConstant(difference);
                                }
                            }

                            phaseDifference_current_index = (phaseDifference_current_index + numSkippedSamples) % phaseDifference.size();
                            break;
                    }
                    last_phase_seqnum = -1;
                    Data_to_display.row(8).head(downsampled_cols - phaseEstParams.edge) = getPhaseDifference_vector();
                }
            }

            if (phasEst_display_all_EEG_channels) {
                for (int i = 0; i < n_EEG_channels_to_use; ++i) {
                    if(i < EEG_channel_names.size()) PhaseEst_channel_names.push_back(EEG_channel_names[i]);
                    else PhaseEst_channel_names.push_back("Undefined");
                }
                emit updatePhaseEstDisplayedData(Data_to_display, triggers_A, triggers_B, triggers_out, downsampled_cols, estimationLength - phaseEstParams.edge);
            } else {
                emit updatePhaseEstDisplayedData(Data_to_display.bottomRows(Data_to_display.rows() - n_EEG_channels_to_use), triggers_A, triggers_B, triggers_out, downsampled_cols, estimationLength - phaseEstParams.edge);
            }
            
            if (spatial_channel_index >= 0 && spatial_channel_index < EEG_channel_names.size()) {
                PhaseEst_channel_names.push_back("Spatial filtered(" + EEG_channel_names[spatial_channel_index] + ")");
            } else {
                PhaseEst_channel_names.push_back("Spatial filtered(Undefined)");
            }

            PhaseEst_channel_names.push_back("Band-pass filtered 9-13Hz & Prediction");
            PhaseEst_channel_names.push_back("Phase angles");
            PhaseEst_channel_names.push_back("Phase difference");
            emit updatePhaseEstwindowNames(PhaseEst_channel_names);

            seq_num_tracker = sequence_number;

            print_debug("Processing end");
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in process_test function: %1").arg(e.what()));
    }
}







// if (performPhaseDifference) {
//                 int numSkippedSamples = (sequence_number - last_phase_seqnum) / downsampling_factor;              // FIGURE OUT A BETTER WAY TO HANDLE DOWNSAMPLING
//                 if (last_phase_seqnum == -1) {
//                     last_phase = phaseAngles(phaseEstParams.edge);
//                     last_phase_seqnum = sequence_number;
//                 } else if (numSkippedSamples > 35) {
//                     phase_diff_hilbert = hilbertTransform(EEG_filter2);
//                     double difference = ang_diff(last_phase, std::arg(phase_diff_hilbert[filter2_length - (numSkippedSamples)]));
                    
//                     // Update the phase difference
//                     if (phaseDifference_current_index + numSkippedSamples < phaseDifference.size()) {
//                         // Safe to modify the segment without wrapping around
//                         phaseDifference.segment(phaseDifference_current_index, numSkippedSamples).setConstant(difference);
//                     } else {
//                         // Calculate the number of elements that fit until the end of the vector
//                         int fitToEnd = phaseDifference.size() - phaseDifference_current_index;

//                         // Update the part that fits
//                         phaseDifference.segment(phaseDifference_current_index, fitToEnd).setConstant(difference);

//                         // Calculate the remaining elements that need to wrap around to the beginning
//                         int overflow = numSkippedSamples - fitToEnd;

//                         // Update the wrapped around part
//                         if (overflow > 0) {
//                             phaseDifference.head(overflow).setConstant(difference);
//                         }
//                     }

//                     phaseDifference_current_index = (phaseDifference_current_index + numSkippedSamples) % phaseDifference.size();
//                     last_phase_seqnum = -1;
//                 }
//                 Data_to_display.row(8).head(downsampled_cols - phaseEstParams.edge) = getPhaseDifference_vector();
//             }