#include "phaseEstimationWorker.h"
#include <iostream>
#include <QThread>


void signalHandlerPhaseEst(int signal) {
    std::cerr << "Error: Segmentation fault (signal " << signal << ")\n";
    std::exit(signal);  // Exit the program
}

phaseEstimationWorker::phaseEstimationWorker(dataHandler &handler, 
                    volatile std::sig_atomic_t &processingWorkerRunning, 
                       phaseEstimateParameters &phaseEstParams_in, 
                                       QObject *parent)
    : QObject(parent), 
      handler(handler), 
      processingWorkerRunning(processingWorkerRunning)
{ 
    setPhaseEstimateParameters(phaseEstParams_in);
}

phaseEstimationWorker::~phaseEstimationWorker()
{ }

void phaseEstimationWorker::setPhaseEstimateParameters(phaseEstimateParameters newParams) {
    std::cout << "New phaseEstimateParameters Parameters: " << newParams << std::endl;
    currentPhaseEstParams = newParams;

    n_channels = handler.get_channel_count();
    samples_to_process = newParams.numberOfSamples;
    downsampling_factor = newParams.downsampling_factor;
    downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;
    EEG_corrected = Eigen::MatrixXd::Zero(n_EEG_channels_to_use, downsampled_cols);
    EEG_spatial = Eigen::VectorXd::Zero(downsampled_cols);

    // Input triggers
    triggers_A = Eigen::VectorXi::Zero(samples_to_process);
    triggers_B = Eigen::VectorXi::Zero(samples_to_process);
    triggers_out = Eigen::VectorXi::Zero(samples_to_process);
    time_stamps = Eigen::VectorXd::Zero(samples_to_process);

    SNR_threshold = newParams.SNR_threshold;
    edge = newParams.edge;
    modelOrder = newParams.modelOrder;
    hilbertWinLength = newParams.hilbertWinLength; 
    stimulation_target = newParams.stimulation_target;
    phase_shift = newParams.phase_shift;
    filter2_length = 250;

    edge_cut_cols = filter2_length - 2 * newParams.edge;
    estimationLength = newParams.edge + std::ceil(newParams.hilbertWinLength / 2);
    display_length = downsampled_cols + estimationLength - newParams.edge;

    EEG_filter2 = Eigen::VectorXd::Zero(filter2_length);
    EEG_predicted.resize(estimationLength, 0.0);
    EEG_hilbert.resize(estimationLength, std::complex<double>(0.0, 0.0));
    phaseAngles = Eigen::VectorXd::Zero(estimationLength);

    phase_diff_hilbert.resize(estimationLength, std::complex<double>(0.0, 0.0));
    phaseDifference = Eigen::VectorXd::Zero(downsampled_cols - newParams.edge);

    outerElectrodeCheckStates_.resize(numOuterElectrodes, true);

    Data_to_display = Eigen::MatrixXd::Zero(9, display_length);
}

void phaseEstimationWorker::handlePreprocessingOutput(const Eigen::MatrixXd &output,
                                           const Eigen::VectorXi &triggers_A_in,
                                           const Eigen::VectorXi &triggers_B_in,
                                           const Eigen::VectorXi &triggers_out_in,
                                           const Eigen::VectorXd &time_stamps_in,
                                           int number_of_samples,
                                           int seq_num) 
{
    EEG_corrected = output;
    triggers_A = triggers_A_in;
    triggers_B = triggers_B_in;
    triggers_out = triggers_out_in;
    time_stamps = time_stamps_in;
    samples_to_process = number_of_samples;
    sequence_number = seq_num;
}

void phaseEstimationWorker::process()
{
    std::signal(SIGSEGV, signalHandlerPhaseEst);

    // Eigen::setNbThreads(std::thread::hardware_concurrency());

    try {
        std::cout << "phaseEstimationWorker start" << '\n';

        // FIR filters
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);

        // Set names for each channel in Data_to_display
        std::vector<std::string> EEG_channel_names;
        std::vector<std::string> PhaseEst_channel_names;
        print_debug("Channel names set");

        int index = 0;
        int last_save_index = -1;
        std::vector<int> trigger_seqNum_list;
        std::vector<int> angle_seqNum_list;
        std::vector<double> angle_list;

        int last_SNR_pass = 0;
        int seq_num_tracker = 0;
        while(processingWorkerRunning) {
            if (processing_pause) {
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
                continue;
            }
            
            print_debug("Processing start");

            // if (sequence_number > 1000000) processingWorkerRunning = false;

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
            
            PhaseEst_channel_names.clear();
            EEG_channel_names = handler.getChannelNames();

            if (EEG_channel_names.size() >= n_EEG_channels_to_use) { 
                std::vector<std::string> EEG_spatial_channel_names(EEG_channel_names.begin(), EEG_channel_names.begin() + n_EEG_channels_to_use);
                emit updateSpatialChannelNames(EEG_spatial_channel_names);
            }

            if (!phaseEstStates.performPhaseEstimation) {
                std::this_thread::sleep_for(std::chrono::milliseconds(5));
                continue;
            }

            emit sendNumSamples(samples_to_process);

            print_debug("Spatial filtering");
            if (spatial_channel_index >= 0 && spatial_channel_index < EEG_corrected.rows()) {
                
                Eigen::VectorXd sum_of_rows = Eigen::VectorXd::Zero(downsampled_cols);
                int outer_channel_index = 0;
                for (int i = 0; i < EEG_corrected.rows(); i++) {
                    if (i == spatial_channel_index) {
                        continue;
                    }
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

            bool SNR_passed = true;
            // SNR check
            if (phaseEstStates.performSNRcheck) {
                print_debug("SNR check");

                double SNR = calculateSNR_max(EEG_spatial.tail(filter2_length), 64, 256, 500.0, 10.0, 6.0);

                if (SNR_max_list.size() < n_SNR_max) {

                    if (!SNR_max_set) {

                        if (SNR_list.size() < n_SNR) {
                            SNR_list.push_back(SNR);
                        } else {
                            SNR_max_temp = *std::max_element(SNR_list.begin(), SNR_list.end());
                            SNR_max_set = true;
                        }

                    } else {

                        SNR_max_list.push_back(SNR_max_temp);

                        //Choose the max or mean of the SNR_max_list for the SNR check
                        // if (SNR_max_final < SNR_max_temp) SNR_max_final = SNR_max_temp;
                        SNR_max_final = std::accumulate(SNR_max_list.begin(), SNR_max_list.end(), 0.0) / SNR_max_list.size();
                        emit sendSNRmax(SNR_max_final);

                        SNR_list.clear();
                        SNR_max_set = false;
                        emit sendSNRmax_list(SNR_max_list);

                        std::this_thread::sleep_for(std::chrono::milliseconds(200));
                    }

                    SNR_passed = false;

                } else {

                    if (SNR < SNR_max_final * SNR_threshold) {
                        SNR_passed = false;
                    }

                }

                // Extra check that ensures the last 200 samples did not pass the SNR check ensuring that we are at the beginning of a clean epoch
                // if (sequence_number < 200 + last_SNR_pass) SNR_passed = false;
            }

            // if (SNR_passed) last_SNR_pass = sequence_number;

            // Filter2
            print_debug("Second filtering");
            // Demean
            EEG_spatial.array() -= EEG_spatial.mean();
            if (phaseEstStates.performFiltering) {
                EEG_filter2 = zeroPhaseLSFIR(EEG_spatial.tail(filter2_length), LSFIR_coeffs_2);
            } else {
                EEG_filter2 = EEG_spatial.tail(filter2_length);
            }

            // Phase estimation
            print_debug("AR predicted");
            if (phaseEstStates.performEstimation) {
                EEG_predicted = fitAndPredictAR_YuleWalker(EEG_filter2.segment(edge, edge_cut_cols), modelOrder, estimationLength);
            }
            
            // Hilbert transform
            print_debug("Hilbert transform");
            if (phaseEstStates.performHilbertTransform) {
                EEG_hilbert = hilbertTransform(EEG_predicted);
            }

            bool save_angle = false;            //Remove later
            int trigger_seqNum = 0;
            // Trigger phase targeting
            print_debug("Phase targeting");
            if (phaseEstStates.performPhaseTargeting) {
                trigger_seqNum = findTargetPhase(EEG_hilbert, phaseAngles, sequence_number, downsampling_factor, edge, phase_shift, stimulation_target);
                if (trigger_seqNum && SNR_passed) { 

                    if (sequence_number > 40000 && trigger_seqNum > 200 + last_save_index) {
                        last_save_index = trigger_seqNum;
                        save_angle = true;
                    }
                    // if (sequence_number > 40000) trigger_seqNum_list.push_back(trigger_seqNum);

                    handler.insertTrigger(trigger_seqNum);
                }
            }

            // Phase histogram
            if (phaseEstStates.performSNRcheck && SNR_max_list.size() < n_SNR_max) {
                int numSkippedSamples = (sequence_number - last_phase_seqnum) / downsampling_factor;
                if (last_phase_seqnum == -1) {
                    last_phase = phaseAngles(edge);
                    last_phase_seqnum = sequence_number;
                } else if (numSkippedSamples > 35) {
                    phase_diff_hilbert = hilbertTransform(EEG_filter2);
                    double difference = ang_diff(last_phase, std::arg(phase_diff_hilbert[filter2_length - numSkippedSamples]));
                    emit polarHistogramAddSample_1(difference);

                    last_phase_seqnum = -1;
                }
            } else if (phaseEstStates.performSNRcheck && PostInitializationCounter < PostInitSamples_n) {
                int numSkippedSamples = (sequence_number - last_phase_seqnum) / downsampling_factor;
                if (last_phase_seqnum == -1 && SNR_passed) {
                    last_phase = phaseAngles(edge);
                    last_phase_seqnum = sequence_number;
                } else if (last_phase_seqnum != -1 && numSkippedSamples > 35) {
                    phase_diff_hilbert = hilbertTransform(EEG_filter2);
                    double difference = ang_diff(last_phase, std::arg(phase_diff_hilbert[filter2_length - numSkippedSamples]));
                    
                    emit polarHistogramAddSample_2(difference);
                    PostInitializationCounter++;
                    last_phase_seqnum = -1;
                }
            }

            Data_to_display.topLeftCorner(5, downsampled_cols) = EEG_corrected;
            
            Data_to_display.row(5).head(downsampled_cols) = EEG_spatial;
            
            Data_to_display.row(6).segment(downsampled_cols - filter2_length, filter2_length - edge) = EEG_filter2.head(filter2_length - edge);


            Data_to_display.row(6).tail(estimationLength) = Eigen::Map<Eigen::VectorXd>(EEG_predicted.data(), EEG_predicted.size());
            
            int phase_length = 32;
            int phase_start = edge - phase_length / 2;
            
            Data_to_display.row(7).segment(downsampled_cols - phase_length / 2, phase_length) = phaseAngles.segment(phase_start, phase_length);

            // Phase difference
            print_debug("Phase difference");
            if (phaseEstStates.performPhaseDifference) {
                int numSkippedSamples = (sequence_number - last_phase_seqnum) / downsampling_factor;
                if (last_phase_seqnum == -1) {
                    last_phase = phaseAngles(edge);
                    last_phase_seqnum = sequence_number;
                } else if (numSkippedSamples > 35) {
                    double difference;
                    switch (phaseErrorType) {
                        case 0:
                            phase_diff_hilbert = hilbertTransform(EEG_filter2);
                            difference = ang_diff(last_phase, std::arg(phase_diff_hilbert[filter2_length - numSkippedSamples]));
                            if (save_angle) {
                                trigger_seqNum_list.push_back(trigger_seqNum);
                                angle_seqNum_list.push_back(sequence_number - 350);
                                angle_list.push_back(std::arg(phase_diff_hilbert[filter2_length - numSkippedSamples]));
                            }

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
                    Data_to_display.row(8).head(downsampled_cols - edge) = getPhaseDifference_vector();
                }
            }

            print_debug("Graph updating");
            if (phaseEstStates.phasEst_display_all_EEG_channels) {
                for (int i = 0; i < n_EEG_channels_to_use; ++i) {
                    if(i < EEG_channel_names.size()) PhaseEst_channel_names.push_back(EEG_channel_names[i]);
                    else PhaseEst_channel_names.push_back("Undefined");
                }
                emit updatePhaseEstDisplayedData(Data_to_display, triggers_A, triggers_B, triggers_out, time_stamps, downsampled_cols, estimationLength - edge);
            } else {
                emit updatePhaseEstDisplayedData(Data_to_display.bottomRows(Data_to_display.rows() - n_EEG_channels_to_use), triggers_A, triggers_B, triggers_out, time_stamps, downsampled_cols, estimationLength - edge);
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
        std::cout << "SNR max: " << SNR_max_final << std::endl;
        std::cout << "Phase estimation finished. Saving data..." << trigger_seqNum_list.size() << std::endl;
        writeMatrixiToCSV("trigger_seqNum_list.csv", vectorToColumnMatrixi(trigger_seqNum_list));
        writeMatrixiToCSV("angle_seqNum_list.csv", vectorToColumnMatrixi(angle_seqNum_list));
        writeMatrixdToCSV("angle_list.csv", vectorToColumnMatrixd(angle_list));
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in processingworker process function: %1").arg(e.what()));
    }
}





