#include "processingworker.h"
#include <iostream>
#include <QThread>


ProcessingWorker::ProcessingWorker(dataHandler &handler, 
                               Eigen::MatrixXd &processed_data, 
                    volatile std::sig_atomic_t &processingWorkerRunning, 
                    const processingParameters &params, 
                                       QObject* parent)
    : QObject(parent), 
      handler(handler), 
      processed_data(processed_data), 
      processingWorkerRunning(processingWorkerRunning), 
      params(params)
{ }

ProcessingWorker::~ProcessingWorker()
{ }

void ProcessingWorker::process()
{
    omp_set_num_threads(10);

    // Eigen::setNbThreads(std::thread::hardware_concurrency());

    try {
        std::cout << "ProcessingWorker start" << '\n';
        std::cout << "params: " << params << '\n';

        // Number of samples to use for processing
        int samples_to_process = params.numberOfSamples;

        //  downsampling
        int downsampling_factor = params.downsampling_factor;
        int downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // removeBCG
        int delay = params.delay;
        int step = (1+2*delay);
        int n_EEG_channels_to_use = 5;
        int n_CWL_channels_to_use = 7;

        // FIR filters
        Eigen::VectorXd LSFIR_coeffs_1;
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_0_80Hz(LSFIR_coeffs_1);
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);

        // phase estimate
        size_t edge = params.edge;
        int edge_cut_cols = downsampled_cols - 2 * edge;
        size_t modelOrder = params.modelOrder;
        size_t hilbertWinLength = params.hilbertWinLength;
        size_t estimationLength = edge + std::ceil(hilbertWinLength / 2);

        // stimulation
        double stimulation_target = params.stimulation_target;
        int phase_shift = params.phase_shift;

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
        Eigen::MatrixXd Data_to_display = Eigen::MatrixXd::Zero(5, downsampled_cols + estimationLength - edge);

        // Set names for each channel in Data_to_display
        std::vector<std::string> processing_channel_names = {"Filter1 & downsample (channel 0)", "removeBCG (channel 0)", "spatial filter", "Filter2 & phase estimation", "phase angle"};
        emit updateProcessingChannelNames(processing_channel_names);

        int seq_num_tracker = 0;
        while(processingWorkerRunning) {
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_process);
            if (seq_num_tracker == sequence_number) continue;

            std::cout << "Samples skipped: " << sequence_number - seq_num_tracker << '\n';
            auto start = std::chrono::high_resolution_clock::now();

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

            // Filter2
            EEG_filter2 = zeroPhaseLSFIR(EEG_spatial, LSFIR_coeffs_2);

            // Phase estimation
            EEG_predicted = fitAndPredictAR_YuleWalker_V2(EEG_filter2.segment(edge, edge_cut_cols), modelOrder, estimationLength);
            
            // Hilbert transform
            EEG_hilbert = hilbertTransform(EEG_predicted);
            
            // Trigger phase targeting
            int trigger_seqNum = findTargetPhase(EEG_hilbert, phaseAngles, sequence_number, downsampling_factor, edge, phase_shift, stimulation_target);
            if (trigger_seqNum) { handler.insertTrigger(trigger_seqNum); }

            seq_num_tracker = sequence_number;
            
            std::chrono::duration<double> total_elapsed = std::chrono::high_resolution_clock::now() - start;
            std::cout << "Time taken: " << total_elapsed.count() << " seconds." << std::endl;


            // Save displayed data and names for each channel
            Data_to_display.row(0).head(downsampled_cols) = EEG_downsampled.row(0);
            Data_to_display.row(1).head(downsampled_cols) = EEG_corrected.row(0);
            Data_to_display.row(2).head(downsampled_cols) = EEG_spatial;
            Data_to_display.row(3).head(downsampled_cols - edge) = EEG_filter2.head(downsampled_cols - edge);

            Data_to_display.row(3).tail(estimationLength) = Eigen::Map<Eigen::VectorXd>(EEG_predicted.data(), EEG_predicted.size());
            Data_to_display.row(4).tail(estimationLength) = phaseAngles.tail(estimationLength);
            processed_data = Data_to_display;
        }
        
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred in process_test function: %1").arg(e.what()));
    }
}











// Testing
void ProcessingWorker::process_timing()
{
    // Set the current thread to use SCHED_RR
    pthread_t this_thread = pthread_self();
    struct sched_param ch_params;
    ch_params.sched_priority = sched_get_priority_max(SCHED_RR);
    
    if (pthread_setschedparam(this_thread, SCHED_RR, &ch_params) != 0) {
        qDebug("Failed to set thread to real-time");
    } else {
        qDebug("Thread set to real-time successfully");
    }

    omp_set_num_threads(5);

    Eigen::setNbThreads(std::thread::hardware_concurrency());

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
        std::vector<double> total_time;

        // Number of samples to use for processing
        int samples_to_process = params.numberOfSamples;

        //  downsampling
        int downsampling_factor = params.downsampling_factor;
        int downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // removeBCG
        int delay = params.delay;
        int step = (1+2*delay);
        int n_CWL_channels_to_use = 7;

        // FIR filters
        Eigen::VectorXd LSFIR_coeffs_1;
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_0_80Hz(LSFIR_coeffs_1);
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);

        // phase estimate
        size_t edge = params.edge;
        int edge_cut_cols = downsampled_cols - 2 * edge;
        size_t modelOrder = params.modelOrder;
        size_t hilbertWinLength = params.hilbertWinLength;
        size_t estimationLength = edge + std::ceil(hilbertWinLength / 2);

        // stimulation
        double stimulation_target = params.stimulation_target;
        int phase_shift = params.phase_shift;

        // Initializing parameters and matrices for algorithms
        int n_channels = handler.get_channel_count();

        Eigen::MatrixXd all_channels = Eigen::MatrixXd::Zero(n_channels, samples_to_process);
        Eigen::MatrixXd EEG_filter1 = Eigen::MatrixXd::Zero(n_channels, samples_to_process);
        Eigen::MatrixXd expCWL = Eigen::MatrixXd::Zero(n_CWL_channels_to_use * step, downsampled_cols);
        Eigen::MatrixXd pinvCWL = Eigen::MatrixXd::Zero(downsampled_cols, downsampled_cols);
        Eigen::MatrixXd EEG_downsampled = Eigen::MatrixXd::Zero(n_channels, downsampled_cols);
        Eigen::MatrixXd EEG_corrected = Eigen::MatrixXd::Zero(n_channels, downsampled_cols);
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
                applyLSFIRFilterMatrix(all_channels * 10, LSFIR_coeffs_1, EEG_filter1);
           
                std::chrono::duration<double> filtering_time = std::chrono::high_resolution_clock::now() - start;
                filtering_total_time += filtering_time.count();

                // Downsampling
                downsample(EEG_filter1, EEG_downsampled, downsampling_factor);

                std::chrono::duration<double> downsampling_time = std::chrono::high_resolution_clock::now() - filtering_time - start;
                downsampling_total_time += downsampling_time.count();

                // CWL
                if (delay > 0) {
                    // Perform delay embedding if delay is positive
                    delayEmbed(EEG_downsampled.middleRows(5, n_CWL_channels_to_use), expCWL, delay);
                } else {
                    expCWL = EEG_downsampled.middleRows(5, n_CWL_channels_to_use);
                }
                removeBCG(EEG_downsampled.middleRows(0, 5), expCWL, pinvCWL, EEG_corrected);
                EEG_spatial = EEG_corrected.row(0) - EEG_corrected.bottomRows(4).colwise().mean();

                std::chrono::duration<double> removeBCG_time = std::chrono::high_resolution_clock::now() - downsampling_time - filtering_time - start;
                removeBCG_total_time += removeBCG_time.count();

                // Second filtering
                EEG_filter2 = zeroPhaseLSFIR(EEG_spatial, LSFIR_coeffs_2);
                // EEG_filter2 = zeroPhaseBW(EEG_spatial, a, b);

                std::chrono::duration<double> pEstFilt_time = std::chrono::high_resolution_clock::now() - removeBCG_time - downsampling_time - filtering_time - start;
                pEstFilt_total_time += pEstFilt_time.count();

                // Phase estimate
                EEG_predicted = fitAndPredictAR_YuleWalker_V2(EEG_filter2.segment(edge, edge_cut_cols), modelOrder, estimationLength);

                std::chrono::duration<double> phaseEstimate_time = std::chrono::high_resolution_clock::now() - pEstFilt_time - removeBCG_time - downsampling_time - filtering_time - start;
                phaseEstimate_total_time += phaseEstimate_time.count();

                // Hilbert transform
                EEG_hilbert = hilbertTransform(EEG_predicted);
                int trigger_seqNum = findTargetPhase(EEG_hilbert, phaseAngles, sequence_number, downsampling_factor, edge, phase_shift, stimulation_target);
                if (trigger_seqNum) { handler.insertTrigger(trigger_seqNum); }

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
                total_time.push_back(total_elapsed.count());
                count++;
            }

            seq_num_tracker = sequence_number;
            QThread::usleep(1);
        }

        if (count > 0) {
            double filtering_average_time = filtering_total_time / count;
            double downsampling_average_time = downsampling_total_time / count;
            double removeBCG_average_time = removeBCG_total_time / count;
            double pEstFilt_average_time = pEstFilt_total_time / count;
            double phaseEstimate_average_time = phaseEstimate_total_time / count;
            double hilbert_average_time = hilbert_total_time / count;
            double average_time = std::reduce(total_time.begin(), total_time.end()) / total_time.size();
            double min_time = * std::min_element(total_time.begin(), total_time.end());
            total_time.erase(total_time.begin());
            double max_time = * std::max_element(total_time.begin(), total_time.end());
            std::cout << "Average filtering time taken: " << filtering_average_time << " seconds." << std::endl;
            std::cout << "Average downsampling time taken: " << downsampling_average_time << " seconds." << std::endl;
            std::cout << "Average removeBCG time taken: " << removeBCG_average_time << " seconds." << std::endl;
            std::cout << "Average second filtering time taken: " << pEstFilt_average_time << " seconds." << std::endl;
            std::cout << "Average phase estimate time taken: " << phaseEstimate_average_time << " seconds." << std::endl;
            std::cout << "Average hilbert time taken: " << hilbert_average_time << " seconds.\n" << std::endl;
            std::cout << "Max total time taken: " << max_time << " seconds." << std::endl;
            std::cout << "Average total time taken: " << average_time << " seconds." << std::endl;
            std::cout << "Min total time taken: " << min_time << " seconds." << std::endl;

            // writeMatrixToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/Betas.csv", Betas);
            writeMatrixdToCSV("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/EEG_corrected.csv", CSV_save);
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
        std::cout << "Pworker start" << '\n';
        std::cout << "params: " << params << '\n';

        // Number of samples to use for processing
        int samples_to_process = params.numberOfSamples;

        //  downsampling
        int downsampling_factor = params.downsampling_factor;
        int downsampled_cols = (samples_to_process + downsampling_factor - 1) / downsampling_factor;

        // FIR filters
        Eigen::VectorXd LSFIR_coeffs_2;
        getLSFIRCoeffs_9_13Hz(LSFIR_coeffs_2);
        int filter2_length = 250;

        // removeBCG
        int delay = params.delay;
        int step = (1+2*delay);
        int n_EEG_channels_to_use = 5;
        int n_CWL_channels_to_use = 7;

        // phase estimate
        size_t edge = params.edge;
        int edge_cut_cols = filter2_length - 2 * edge;
        size_t modelOrder = params.modelOrder;
        size_t hilbertWinLength = params.hilbertWinLength;
        size_t estimationLength = edge + std::ceil(hilbertWinLength / 2);

        // stimulation
        double stimulation_target = params.stimulation_target;
        int phase_shift = params.phase_shift;

        // Initializing parameters and matrices for algorithms
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
        Eigen::MatrixXd Data_to_display = Eigen::MatrixXd::Zero(5, downsampled_cols + estimationLength - edge);

        Eigen::MatrixXd CSV_CWL_save = Eigen::MatrixXd::Zero(164, downsampled_cols);
        Eigen::MatrixXd CSV_filter2_save = Eigen::MatrixXd::Zero(164, filter2_length);
        Eigen::MatrixXd CSV_predict_save = Eigen::MatrixXd::Zero(164, estimationLength);
        int csv_save_count = 0;
        int last_processed_index = 0;
        std::vector<int> index_list;
        
        Eigen::VectorXd EEG_filter2_total = Eigen::VectorXd::Zero(downsampled_cols);
        std::vector<std::complex<double>> EEG_hilbert_total(estimationLength, std::complex<double>(0.0, 0.0));
        Eigen::VectorXd phaseAngles_total = Eigen::VectorXd::Zero(estimationLength);
        Eigen::MatrixXd phaseAngles_out = Eigen::MatrixXd::Zero(164, estimationLength);
        Eigen::MatrixXd phaseAngles_difference = Eigen::MatrixXd::Zero(164, estimationLength);

        // Set names for each channel in Data_to_display
        std::vector<std::string> processing_channel_names = {"Filter1 & downsample (channel 0)", "removeBCG (channel 0)", "spatial filter", "Filter2 & phase estimation", "phase angle"};
        emit updateProcessingChannelNames(processing_channel_names);
        
        int seq_num_tracker = 0;
        while(processingWorkerRunning) {
            int sequence_number = handler.getLatestDataInOrder(all_channels, samples_to_process);
            if (/*seq_num_tracker == sequence_number ||*/last_processed_index == sequence_number || !(sequence_number % 10000 == 0)) {
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

            // Filter2
            EEG_filter2 = zeroPhaseLSFIR(EEG_spatial.segment(EEG_spatial.size() / 2 - filter2_length, filter2_length), LSFIR_coeffs_2);
            EEG_filter2_total = zeroPhaseLSFIR(EEG_spatial, LSFIR_coeffs_2);
            CSV_filter2_save.row(csv_save_count) = EEG_filter2;

            // Phase estimation
            EEG_predicted = fitAndPredictAR_YuleWalker_V2(EEG_filter2.segment(edge, edge_cut_cols), modelOrder, estimationLength);

            // Hilbert transform
            EEG_hilbert = hilbertTransform(EEG_predicted);
            Eigen::VectorXd EEG_filter2_cut = EEG_filter2_total.segment(EEG_spatial.size() / 2 - edge, estimationLength);
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
            index_list.push_back(trigger_seqNum);
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
            processed_data = Data_to_display;
            csv_save_count++;
            last_processed_index = sequence_number;
        }

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







