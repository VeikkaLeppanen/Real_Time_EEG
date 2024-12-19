#include "dataHandler.h"

// Resets the buffers for new channelcount and sampling rates
void dataHandler::reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate) {

    channel_count_ = channel_count;
    sampling_rate_ = sampling_rate;
    simulation_delivery_rate_ = simulation_delivery_rate;

    buffer_capacity_ = buffer_length_in_seconds_ * sampling_rate;
    current_data_index_ = 0;

    {
        std::lock_guard<std::mutex> lock(this->dataMutex);
        sample_buffer_ = Eigen::MatrixXd::Zero(channel_count, buffer_capacity_);
        time_stamp_buffer_ = Eigen::VectorXd::Zero(buffer_capacity_);
        trigger_buffer_A = Eigen::VectorXi::Zero(buffer_capacity_);
        trigger_buffer_B = Eigen::VectorXi::Zero(buffer_capacity_);
        trigger_buffer_out = Eigen::VectorXi::Zero(buffer_capacity_);
    }

    processing_sample_vector = Eigen::VectorXd::Zero(channel_count);

    GACorr_ = GACorrection(channel_count, GA_average_length, TA_length);

    RTfilter_.reset_filter(channel_count);

    handler_state = WAITING_FOR_STOP;
}
    
void dataHandler::reset_GACorr(int TA_length_input, int GA_average_length_input) {
    TA_length = TA_length_input; 
    GA_average_length = GA_average_length_input;
    GACorr_ = GACorrection(channel_count_, GA_average_length, TA_length);
    TA_tracker = 10000000;
}
















// Buffer functions
// Add a single sample to each channel
void dataHandler::addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger_A, const int &trigger_B, const int &SeqNo) {

    try {
        // Debug: Validate sample size
        if (samples.size() == 0) {
            std::cerr << "Error: Empty samples vector." << std::endl;
            return;
        }
        std::lock_guard<std::mutex> lock(this->dataMutex);

        if (trigger_A == 1) { 
            TA_tracker = SeqNo;
        }

        // Gradient artifact correction
        if (Apply_GACorr && TA_tracker > 0 && SeqNo >= TA_tracker && SeqNo < TA_tracker + TA_length) {

            int temp_index = SeqNo - TA_tracker;

            if (temp_index < 0) {
                std::cerr << "Error: temp_index is negative." << std::endl;
                return;
            }

            if (temp_index >= GACorr_.getTemplateSize()) {
                std::cerr << "Error: temp_index exceeds GACorr template size." << std::endl;
                return;
            }

            processing_sample_vector = samples - GACorr_.getTemplateCol(temp_index);
            
            GACorr_.update_template(temp_index, samples);
        } else {
            processing_sample_vector = samples;
        }

        // Baseline average
        if (Apply_baseline) {
            if (baseline_length > 0) {
                // Ensure baseline_average and baseline_matrix are initialized correctly
                if (baseline_average.size() != samples.size()) {
                    std::cerr << "Warning: Resetting baseline structures due to size mismatch." << std::endl;
                    baseline_average = Eigen::VectorXd::Zero(samples.size());
                    baseline_matrix = Eigen::MatrixXd::Zero(samples.size(), baseline_length);
                    baseline_index = 0;
                }

                // Update baseline_average
                baseline_average += (processing_sample_vector - baseline_matrix.col(baseline_index)) / baseline_length;
                baseline_matrix.col(baseline_index) = processing_sample_vector;

                // Avoid dividing by zero in baseline_average
                for (int i = 0; i < baseline_average.size(); ++i) {
                    if (baseline_average(i) == 0) {
                        baseline_average(i) = 1; // Prevent division by zero
                    }
                }

                // Update processing_sample_vector
                processing_sample_vector = processing_sample_vector.cwiseQuotient(baseline_average);

                // Update baseline_index
                baseline_index = (baseline_index + 1) % baseline_length;
            } else {
                std::cerr << "Error: baseline_length is zero or negative." << std::endl;
            }
        }

        // Filtering
        if (Apply_filter) {
            processing_sample_vector = RTfilter_.processSample(processing_sample_vector);
        }

        // Debug: Check buffer capacity before indexing
        if (current_data_index_ >= buffer_capacity_) {
            std::cerr << "Error: current_data_index exceeds buffer capacity." << std::endl;
            return;
        }

        // Triggering
        if (getTriggerEnableStatus() && shouldTrigger(SeqNo) && checkTimeLimit()) {
            latest_trigger_time = std::chrono::system_clock::now();
            
            if (getTriggerConnectStatus()) {
                switch (TMS_connectionType) {
                    case COM:
                        send_trigger();
                        break;
                    case TTL:
                        send_trigger_TTL();
                        break;
                }
            }
            
            seqNum_list.push_back(SeqNo);

            removeTrigger(SeqNo);
            trigger_buffer_out(current_data_index_) = 1;
        } else {
            trigger_buffer_out(current_data_index_) = 0;
        }

        // if (SeqNo >= 1600000 && !data_saved) {
        //     std::cout << "Saving data..." << seqNum_list.size() << std::endl;
        //     writeMatrixiToCSV("trigger_seqNum_list.csv", vectorToColumnMatrixi(seqNum_list));
        //     data_saved = true;
        // }

        
        // Samples, timestamp, triggers, and buffer index update
        sample_buffer_.col(current_data_index_) = processing_sample_vector;
        time_stamp_buffer_(current_data_index_) = time_stamp;
        trigger_buffer_A(current_data_index_) = trigger_A;
        trigger_buffer_B(current_data_index_) = trigger_B;
        current_sequence_number_ = SeqNo;
        current_data_index_ = (current_data_index_ + 1) % buffer_capacity_;

        new_data_available = true;
        

        // Notify the processing thread
        data_condition.notify_one();

        // if (SAVE_INDEX_TRACKER >= 1000000 && !data_saved) {
        //     data_saved = true;
        //     writeMatrixdToCSV("data_GAcorr_filter1_interleaved.csv", save_matrix);
        // }
        // if (SAVE_INDEX_TRACKER <= 1000000) {
        //     save_matrix.col(SAVE_INDEX_TRACKER) = processing_sample_vector;
        //     SAVE_INDEX_TRACKER++;
        // }

    } catch (const std::exception& e) {
        std::cerr << "Datahandler exception: " << e.what() << '\n';
        std::cerr << boost::stacktrace::stacktrace();  // Ensure boost stacktrace is linked correctly
    }
}

int dataHandler::getLatestDataAndTriggers(Eigen::MatrixXd &output, 
                                          Eigen::VectorXi &triggers_A, 
                                          Eigen::VectorXi &triggers_B, 
                                          Eigen::VectorXi &triggers_out, 
                                          Eigen::VectorXd &time_stamps, 
                                                      int number_of_samples) {

    if (last_retrieved_sequence_number_ == current_sequence_number_) {
        return current_sequence_number_;
    }

    // Ensure matrices are resized correctly
    if (output.rows() != channel_count_ || output.cols() != number_of_samples) output.resize(channel_count_, number_of_samples);
    if (triggers_A.cols() != number_of_samples) triggers_A.resize(number_of_samples);
    if (triggers_B.cols() != number_of_samples) triggers_B.resize(number_of_samples);
    if (triggers_out.cols() != number_of_samples) triggers_out.resize(number_of_samples);
    if (time_stamps.cols() != number_of_samples) time_stamps.resize(number_of_samples);

    // Calculate the number of samples that fit before reaching the end of the buffer
    int fitToEnd = std::min(number_of_samples, static_cast<int>(current_data_index_));
    int overflow = number_of_samples - fitToEnd;

    std::lock_guard<std::mutex> lock(this->dataMutex); // Protect shared data access

    // Debugging check before accessing rightCols and middleCols
    if (fitToEnd > 0) {
        output.rightCols(fitToEnd) = sample_buffer_.middleCols(current_data_index_ - fitToEnd, fitToEnd);
        triggers_A.tail(fitToEnd) = trigger_buffer_A.segment(current_data_index_ - fitToEnd, fitToEnd);
        triggers_B.tail(fitToEnd) = trigger_buffer_B.segment(current_data_index_ - fitToEnd, fitToEnd);
        triggers_out.tail(fitToEnd) = trigger_buffer_out.segment(current_data_index_ - fitToEnd, fitToEnd);
        time_stamps.tail(fitToEnd) = time_stamp_buffer_.segment(current_data_index_ - fitToEnd, fitToEnd);
    }

    // Debugging check before accessing leftCols and middleCols
    if (overflow > 0) {
        output.leftCols(overflow) = sample_buffer_.middleCols(buffer_capacity_ - overflow, overflow);
        triggers_A.head(overflow) = trigger_buffer_A.segment(buffer_capacity_ - overflow, overflow);
        triggers_B.head(overflow) = trigger_buffer_B.segment(buffer_capacity_ - overflow, overflow);
        triggers_out.head(overflow) = trigger_buffer_out.segment(buffer_capacity_ - overflow, overflow);
        time_stamps.head(overflow) = time_stamp_buffer_.segment(buffer_capacity_ - overflow, overflow);
    }

    last_retrieved_sequence_number_ = current_sequence_number_;
    return current_sequence_number_;
}







// MAGPRO FUNCTIONS

int dataHandler::connectTriggerPort() {
    return magPro_3G.connectTriggerPort();
}

void dataHandler::send_trigger() {
    magPro_3G.trig();
}

void dataHandler::set_enable(bool status) {
    magPro_3G.set_enable(status);
    triggerEnableState = status;
}

void dataHandler::set_amplitude(int amplitude) {
    magPro_3G.set_amplitude(amplitude);
}

void dataHandler::magPro_set_mode(int mode, int direction, int waveform, int burst_pulses, float ipi, float ba_ratio, bool delay) {
    magPro_3G.set_mode(mode, direction, waveform, burst_pulses, ipi, ba_ratio, delay);
}

void dataHandler::magPro_request_mode_info() {
    magPro_3G.request_G3_to_send_mode_info();
    magPro_3G.sleep(0.3);
    magPro_3G.handle_input_queue();
    magPro_3G.sleep(0.3);
}

void dataHandler::get_mode_info(int &mode, int &direction, int &waveform, int &burst_pulses, float &ipi, float &ba_ratio, bool &enabled) {
    mode = magPro_3G.get_current_mode();
    direction = magPro_3G.get_current_direction();
    waveform = magPro_3G.get_current_waveform();
    burst_pulses = magPro_3G.get_current_burst_pulses();
    ipi = magPro_3G.get_current_ipi();
    ba_ratio = magPro_3G.get_current_ba_ratio();
    enabled = magPro_3G.get_current_enabled();
}

// TTL functions

int dataHandler::connectTriggerPort_TTL() {
    if (labjack_handle == -1) {
        int err = LJM_Open(LJM_dtANY, LJM_ctANY, "LJM_idANY", &labjack_handle);
        if (err != LJME_NOERROR) {
            labjack_handle = -1;
            std::cerr << "Failed to connect to LabJack" << std::endl;
        } else {
            // PrintDeviceInfoFromHandle(labjack_handle);
            std::cerr << "Successfully connected to LabJack." << std::endl;
            return 0;
        }
    }
    return 1;
}

void dataHandler::send_trigger_TTL() {
  if (labjack_handle == -1) {
    std::cerr << "LabJack is not connected. Skipping trigger.";
    return;
  }

  /* Set output port state to high. */
  int err = LJM_eWriteName(labjack_handle, "FIO4", 1);
  if (err != LJME_NOERROR) {
    std::cerr << "LabJack failed to set output port high." << std::endl;
    return;
  }

  /* Wait for one millisecond. */
//   std::this_thread::sleep_for(std::chrono::milliseconds(1));

  /* Set output port state to low. */
  err = LJM_eWriteName(labjack_handle, "FIO4", 0);
  if (err != LJME_NOERROR) {
    std::cerr << "LabJack failed to set output port low." << std::endl;
    return;
  }

}

void dataHandler::set_enable_TTL(bool status) {
    triggerEnableState = status;
}
