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

            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_process);
            
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
            print_debug("Delay Embedding");
            if (delay > 0) {
                // Perform delay embedding if delay is positive
                delayEmbed(EEG_downsampled.middleRows(5, n_CWL_channels_to_use), expCWL, delay);
            } else {
                expCWL = EEG_downsampled.middleRows(5, n_CWL_channels_to_use);
            }

            print_debug("removeBCG");
            removeBCG(EEG_downsampled.middleRows(0, 5), expCWL, pinvCWL, EEG_corrected);
            
            // Update the EEG window graph data
            switch (display_state) 
            {
                case displayState::RAW:
                    emit updateEEGDisplayedData(all_channels);
                    emit updateEEGwindowNames(EEG_channel_names);
                    break;
                case displayState::DOWNSAMPLED:
                    emit updateEEGDisplayedData(EEG_downsampled);
                    emit updateEEGwindowNames(EEG_channel_names);
                    break;
                case displayState::CWL:
                    emit updateEEGDisplayedData(EEG_corrected);
                    emit updateEEGwindowNames(EEG_channel_names);
                    break;
            }

            if (!performPhaseEstimation) continue;

            print_debug("Spatial filtering");
            if (spatial_channel_index >= 0 && spatial_channel_index < EEG_corrected.rows()) {
                
                Eigen::VectorXd sum_of_all_rows = EEG_corrected.colwise().sum();
                sum_of_all_rows -= EEG_corrected.row(spatial_channel_index).transpose();

                Eigen::VectorXd mean_of_remaining = sum_of_all_rows / (EEG_corrected.rows() - 1);
                EEG_spatial = EEG_corrected.row(spatial_channel_index) - mean_of_remaining.transpose();

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
            Data_to_display.row(7).tail(estimationLength) = phaseAngles.tail(estimationLength);

            if (phasEst_display_all_EEG_channels) {
                for (int i = 0; i < n_EEG_channels_to_use; ++i) {
                    if(i < EEG_channel_names.size()) PhaseEst_channel_names.push_back(EEG_channel_names[i]);
                    else PhaseEst_channel_names.push_back("Undefined");
                }
                emit updatePhaseEstDisplayedData(Data_to_display);
            } else {
                emit updatePhaseEstDisplayedData(Data_to_display.bottomRows(Data_to_display.rows() - n_EEG_channels_to_use));
            }
            
            if (spatial_channel_index >= 0 && spatial_channel_index < EEG_channel_names.size()) {
                PhaseEst_channel_names.push_back("Spatial filtered(" + EEG_channel_names[spatial_channel_index] + ")");
            } else {
                PhaseEst_channel_names.push_back("Spatial filtered(Undefined)");
            }

            PhaseEst_channel_names.push_back("Band-pass filtered 9-13Hz & Prediction");
            PhaseEst_channel_names.push_back("Phase angles");
            emit updatePhaseEstwindowNames(PhaseEst_channel_names);

            seq_num_tracker = sequence_number;

            print_debug("Processing end");
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in process_test function: %1").arg(e.what()));
    }
}







void ProcessingWorker::process_testing()
{
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
        int downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // removeBCG
        int delay = prepParams.delay;
        int step = (1+2*delay);
        int n_EEG_channels_to_use = 5;
        int n_CWL_channels_to_use = 7;

        // FIR filters
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);
        int filter2_length = 500;

        // phase estimate
        size_t edge = phaseEstParams.edge;
        int edge_cut_cols = downsampled_cols - 2 * edge;
        int predict_start_index = downsampled_cols / 2;
        size_t modelOrder = phaseEstParams.modelOrder;
        size_t hilbertWinLength = phaseEstParams.hilbertWinLength;
        size_t estimationLength = edge + std::ceil(hilbertWinLength / 2);

        // stimulation
        double stimulation_target = phaseEstParams.stimulation_target;
        int phase_shift = phaseEstParams.phase_shift;
        int stim_window_delay = 2500;

        int n_channels = handler.get_channel_count();

        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_channels, samples_to_process);
        Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_channels, downsampled_cols);
        Eigen::MatrixXd expCWL = Eigen::MatrixXd::Zero(n_CWL_channels_to_use * step, downsampled_cols);
        Eigen::MatrixXd pinvCWL = Eigen::MatrixXd::Zero(downsampled_cols, downsampled_cols);
        Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_EEG_channels_to_use, downsampled_cols);
        Eigen::VectorXd EEG_spatial = Eigen::VectorXd::Zero(downsampled_cols);
        Eigen::VectorXd EEG_filter2 = Eigen::VectorXd::Zero(filter2_length);
        std::vector<double> EEG_predicted(estimationLength, 0.0);
        std::vector<std::complex<double>> EEG_hilbert(estimationLength, std::complex<double>(0.0, 0.0));
        Eigen::VectorXd phaseAngles = Eigen::VectorXd::Zero(estimationLength);

        Eigen::VectorXd EEG_predicted_EIGEN = Eigen::VectorXd::Zero(EEG_predicted.size());
        Eigen::MatrixXd Data_to_display = Eigen::MatrixXd::Zero(6, downsampled_cols + estimationLength - edge);

        int number_of_save_rows = 1000;
        Eigen::MatrixXd CSV_CWL_save = Eigen::MatrixXd::Zero(number_of_save_rows, downsampled_cols);
        Eigen::MatrixXd CSV_filter2_save = Eigen::MatrixXd::Zero(number_of_save_rows, filter2_length);
        Eigen::MatrixXd CSV_predict_save = Eigen::MatrixXd::Zero(number_of_save_rows, estimationLength);
        int csv_save_count = 0;
        int last_processed_index = 0;
        std::vector<int> index_list;

        double total_SNR = 0;
        int snr_count = 0;

        Eigen::VectorXd EEG_filter2_total = Eigen::VectorXd::Zero(downsampled_cols);
        std::vector<std::complex<double>> EEG_hilbert_total(estimationLength, std::complex<double>(0.0, 0.0));
        Eigen::VectorXd phaseAngles_total = Eigen::VectorXd::Zero(estimationLength);
        Eigen::MatrixXd phaseAngles_out = Eigen::MatrixXd::Zero(number_of_save_rows, estimationLength);
        Eigen::MatrixXd phaseAngles_difference = Eigen::MatrixXd::Zero(number_of_save_rows, estimationLength);

        // Set names for each channel in Data_to_display
        std::vector<std::string> processing_channel_names = {"Filter1 & downsample (channel 0)", "removeBCG (channel 0)", "spatial filter", "Filter2 & phase estimation", "phase angle"};
        emit updatePhaseEstwindowNames(processing_channel_names);
        
        int seq_num_tracker = 0;
        while(processingWorkerRunning) {
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_process);
            if (/*seq_num_tracker == sequence_number ||*/last_processed_index == sequence_number || !(sequence_number % 1000 == 0) || sequence_number < 10001) {
                std::this_thread::sleep_for(std::chrono::microseconds(1));
                continue;
            }

            std::cout << "Samples skipped: " << sequence_number - seq_num_tracker << '\n';
            auto start = std::chrono::high_resolution_clock::now();

            // Downsampling
            downsample(all_channels, EEG_downsampled, downsampling_factor);

            // CWL
            if (delay > 0) {
                delayEmbed(EEG_downsampled.middleRows(5, n_CWL_channels_to_use), expCWL, delay);
            } else {
                expCWL = EEG_downsampled.middleRows(5, n_CWL_channels_to_use);
            }
            removeBCG(EEG_downsampled.middleRows(0, 5), expCWL, pinvCWL, EEG_corrected);
            CSV_CWL_save.row(csv_save_count) = EEG_corrected.row(0);
            EEG_spatial = EEG_corrected.row(0) - EEG_corrected.bottomRows(4).colwise().mean();

            // SNR check
            double SNR = calculateSNR(EEG_spatial.segment(predict_start_index - filter2_length, filter2_length), 64, 500, 500.0, 10.0, 2.0);
            total_SNR += SNR;
            snr_count++;
            if (SNR < 10.0) {
                std::cout << "SNR too small: " << SNR << '\n';
                continue;
            }

            // Filter2
            EEG_filter2 = zeroPhaseLSFIR(EEG_spatial.segment(predict_start_index - filter2_length, filter2_length), LSFIR_coeffs_2);
            EEG_filter2_total = zeroPhaseLSFIR(EEG_spatial, LSFIR_coeffs_2);
            CSV_filter2_save.row(csv_save_count) = EEG_filter2;

            // Phase estimation
            EEG_predicted = fitAndPredictAR_YuleWalker_V2(EEG_filter2.segment(edge, edge_cut_cols), modelOrder, estimationLength);

            // Hilbert transform
            EEG_hilbert = hilbertTransform(EEG_predicted);
            Eigen::VectorXd EEG_filter2_cut = EEG_filter2_total.segment(predict_start_index - edge, estimationLength);
            std::vector<double> dtf(EEG_filter2_cut.data(), EEG_filter2_cut.data() + EEG_filter2_cut.size());
            EEG_hilbert_total = hilbertTransform(dtf);

            for (std::size_t i = 0; i < EEG_hilbert.size(); ++i) {
                phaseAngles(i) = std::arg(EEG_hilbert[i]);
                phaseAngles_out(csv_save_count, i) = std::arg(EEG_hilbert_total[i]);
                phaseAngles_total(i) = std::arg(EEG_hilbert_total[i]);
                // phaseAngles_total(i) = stimulation_target;
            }

            phaseAngles_difference.row(csv_save_count) = ang_diff(phaseAngles_total, phaseAngles);

            // Trigger phase targeting
            int trigger_seqNum = findTargetPhase(EEG_hilbert, phaseAngles, sequence_number, downsampling_factor, edge, phase_shift, stimulation_target);
            index_list.push_back(trigger_seqNum - sequence_number);
            // if (trigger_seqNum) { 
            //     handler.insertTrigger(trigger_seqNum);
            //     // if (sequence_number % 1000 == 0) { index_list.push_back(trigger_seqNum); }
            // }

            CSV_predict_save.row(csv_save_count) = Eigen::Map<Eigen::VectorXd>(EEG_predicted.data(), EEG_predicted.size());

            seq_num_tracker = sequence_number;

            std::chrono::duration<double> total_elapsed = std::chrono::high_resolution_clock::now() - start;
            std::cout << "Time taken: " << total_elapsed.count() << " sec." << " Seq_num: " << sequence_number << " , count: " << csv_save_count << std::endl;


            // Save displayed data
            Data_to_display.row(0).head(downsampled_cols) = EEG_downsampled.row(0);
            Data_to_display.row(1).head(downsampled_cols) = EEG_corrected.row(0);
            Data_to_display.row(2).head(downsampled_cols) = EEG_spatial;
            Data_to_display.row(3).segment(downsampled_cols - filter2_length, filter2_length - edge) = EEG_filter2.head(filter2_length - edge);

            // EEG_predicted_EIGEN = Eigen::Map<Eigen::VectorXd>(EEG_predicted.data(), EEG_predicted.size());
            Data_to_display.row(3).tail(estimationLength) = Eigen::Map<Eigen::VectorXd>(EEG_predicted.data(), EEG_predicted.size());
            Data_to_display.row(4).tail(estimationLength) = phaseAngles.tail(estimationLength);

            {
                std::lock_guard<std::mutex> lock(this->dataMutex); // Protect shared data access
                processed_data = Data_to_display;
            }

            csv_save_count++;
            last_processed_index = sequence_number;
            std::cout << "End of processing loop " << std::endl;
        }
        
        std::cout << "Total SNR: " << total_SNR / snr_count << '\n';
        writeMatrixdToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/CWL.csv", CSV_CWL_save);
        writeMatrixdToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/filter2_end.csv", CSV_filter2_save);
        writeMatrixdToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/predicted_end.csv", CSV_predict_save);
        writeMatrixdToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/phaseAngles_total.csv", phaseAngles_out);
        writeMatrixdToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/phaseAngles_difference.csv", phaseAngles_difference);
        writeMatrixiToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/index_list.csv", vectorToColumnMatrixi(index_list));
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in process_test function: %1").arg(e.what()));
        std::cerr << "Processing exception: " << e.what() << '\n';
        std::cerr << boost::stacktrace::stacktrace();
    }
}







void ProcessingWorker::process_stimulation_testing()
{
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
        int downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // removeBCG
        int delay = prepParams.delay;
        int step = (1+2*delay);
        int n_EEG_channels_to_use = 5;
        int n_CWL_channels_to_use = 7;

        // FIR filters
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);

        // phase estimate
        size_t edge = phaseEstParams.edge;
        int edge_cut_cols = downsampled_cols - 2 * edge;
        size_t modelOrder = phaseEstParams.modelOrder;
        size_t hilbertWinLength = phaseEstParams.hilbertWinLength;
        size_t estimationLength = edge + std::ceil(hilbertWinLength / 2);

        // stimulation
        double stimulation_target = phaseEstParams.stimulation_target;
        int phase_shift = phaseEstParams.phase_shift;
        int stim_window_delay = 2500;

        int n_channels = handler.get_channel_count();

        // Memory preallocation for processing matrices
        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_channels, samples_to_process);
        Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_channels, downsampled_cols);
        Eigen::MatrixXd expCWL = Eigen::MatrixXd::Zero(n_CWL_channels_to_use * step, downsampled_cols);
        Eigen::MatrixXd pinvCWL = Eigen::MatrixXd::Zero(downsampled_cols, downsampled_cols);
        Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_EEG_channels_to_use, downsampled_cols);
        Eigen::VectorXd EEG_spatial = Eigen::VectorXd::Zero(downsampled_cols);
        Eigen::VectorXd EEG_filter2 = Eigen::VectorXd::Zero(downsampled_cols);
        std::vector<double> EEG_predicted(estimationLength, 0.0);
        std::vector<std::complex<double>> EEG_hilbert(estimationLength, std::complex<double>(0.0, 0.0));
        Eigen::VectorXd phaseAngles(estimationLength);

        // Data visualization matrices
        Eigen::VectorXd EEG_predicted_EIGEN = Eigen::VectorXd::Zero(EEG_predicted.size());
        int display_length = downsampled_cols + estimationLength - edge;
        Eigen::MatrixXd Data_to_display = Eigen::MatrixXd::Zero(6, display_length);

        // Set names for each channel in Data_to_display
        std::vector<std::string> processing_channel_names = {"Filter1 & downsample (channel 0)", "removeBCG (channel 0)", "spatial filter", "Filter2 & phase estimation", "phase angle", "Stimulation"};
        emit updatePhaseEstwindowNames(processing_channel_names);

        int seq_num_tracker = 0;
        int stimulation_tracker = -1;
        while(processingWorkerRunning) {
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_process);
            if (seq_num_tracker == sequence_number) continue;

            // std::cout << "Samples skipped: " << sequence_number - seq_num_tracker << '\n';
            // auto start = std::chrono::high_resolution_clock::now();

            // Downsampling
            downsample(all_channels, EEG_downsampled, downsampling_factor);

            // CWL
            if (delay > 0) {
                // Perform delay embedding if delay is positive
                delayEmbed(EEG_downsampled.middleRows(5, n_CWL_channels_to_use), expCWL, delay);
            } else {
                expCWL = EEG_downsampled.middleRows(5, n_CWL_channels_to_use);
            }
            removeBCG(EEG_downsampled.middleRows(0, 5), expCWL, pinvCWL, EEG_corrected);
            EEG_spatial = EEG_corrected.row(0) - EEG_corrected.bottomRows(4).colwise().mean();

            // SNR check
            double SNR = calculateSNR(EEG_spatial, 64, 500, 500.0, 10.0, 2.0);
            if (SNR < 2) {
                std::cout << "SNR too small: " << SNR << '\n';
                continue;
            }

            // Filter2
            EEG_filter2 = zeroPhaseLSFIR(EEG_spatial, LSFIR_coeffs_2);

            // Phase estimation
            EEG_predicted = fitAndPredictAR_YuleWalker_V2(EEG_filter2.segment(edge, edge_cut_cols), modelOrder, estimationLength);
            
            // Hilbert transform
            EEG_hilbert = hilbertTransform(EEG_predicted);
            
            // Trigger phase targeting
            int trigger_seqNum = findTargetPhase(EEG_hilbert, phaseAngles, sequence_number, downsampling_factor, edge, phase_shift, stimulation_target);
            if (trigger_seqNum) { 
                handler.insertTrigger(trigger_seqNum); 
                if (stimulation_tracker < 0) stimulation_tracker++;
            }

            // Trigger visualization
            if (stimulation_tracker > stim_window_delay) {
                int stim_length = (stim_window_delay + stimulation_tracker) / downsampling_factor;
                int stim_length_cut = (stim_window_delay * 2) / downsampling_factor;
                Data_to_display.row(5).segment(downsampled_cols - stim_length_cut, stim_length_cut) = EEG_filter2.tail(stim_length).head(stim_length_cut);
                stimulation_tracker = -1;
            } else if (stimulation_tracker > -1) {
                stimulation_tracker += sequence_number - seq_num_tracker;
            }

            seq_num_tracker = sequence_number;
            
            // std::chrono::duration<double> total_elapsed = std::chrono::high_resolution_clock::now() - start;
            // std::cout << "Time taken: " << total_elapsed.count() << " seconds." << std::endl;


            // Save displayed data and names for each channel
            Data_to_display.row(0).head(downsampled_cols) = EEG_downsampled.row(0);
            Data_to_display.row(1).head(downsampled_cols) = EEG_corrected.row(0);
            Data_to_display.row(2).head(downsampled_cols) = EEG_spatial;
            Data_to_display.row(3).head(downsampled_cols - edge) = EEG_filter2.head(downsampled_cols - edge);

            Data_to_display.row(3).tail(estimationLength) = Eigen::Map<Eigen::VectorXd>(EEG_predicted.data(), EEG_predicted.size());
            Data_to_display.row(4).tail(estimationLength) = phaseAngles.tail(estimationLength);
            {
                std::lock_guard<std::mutex> lock(this->dataMutex); // Protect shared data access
                processed_data = Data_to_display;
            }
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in process_test function: %1").arg(e.what()));
    }
}