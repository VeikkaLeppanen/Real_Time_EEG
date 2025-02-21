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
        
        PreProcessing_output_save = Eigen::MatrixXd::Zero(number_of_EEG_channels, buffer_capacity_);
        for (int i = 0; i < channel_count; i++) {
            channel_names_.push_back("unknown_channel_" + std::to_string(i + 1));
        }
    }

    sample_buffer_save = Eigen::MatrixXd::Zero(20, buffer_capacity_);

    processing_sample_vector = Eigen::VectorXd::Zero(channel_count);

    baseline_average = Eigen::VectorXd::Zero(channel_count);
    GA_sum = Eigen::VectorXd::Zero(channel_count);
    GA_average = Eigen::VectorXd::Zero(channel_count);
    baseline_correction = Eigen::VectorXd::Zero(channel_count);

    GACorr_ = GACorrection(channel_count, GA_average_length, GA_shift_front + TA_length + GA_shift_back);

    RTfilter_.reset_filter(channel_count);

    handler_state = WAITING_FOR_STOP;
}
    
void dataHandler::reset_GACorr(int TA_length_input, int GA_average_length_input) {
    TA_length = TA_length_input; 
    GA_average_length = GA_average_length_input;
    GACorr_ = GACorrection(channel_count_, GA_average_length, GA_shift_front + TA_length + GA_shift_back);
    TA_tracker = -1;
    GA_tracker = -1;
    if (TA_length == TR_length) {
        GA_is_continuous = true;
        GA_shift_front = 0;
        GA_shift_back = 0;
    } else {
        GA_is_continuous = false;
        GA_shift_front = 100;
        GA_shift_back = 100;
    }
}
















// Buffer functions
/*
Add a single sample to each channel. This fuction also handles some preprocessing such as TA, TR, GA, baseline correction, geometric sum correction, and filtering.
Realtime operations such as triggering are also handled here.
*/
void dataHandler::addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger_A, const int &trigger_B, const int &SeqNo) {
    auto start = std::chrono::high_resolution_clock::now();

    try {
        // Debug: Validate sample size
        if (samples.size() == 0) {
            std::cerr << "Error: Empty samples vector." << std::endl;
            return;
        }
        std::lock_guard<std::mutex> lock(this->dataMutex);

        if (trigger_A == 1) { 
            TA_tracker = SeqNo;
            if (GA_is_continuous) GA_tracker = SeqNo;
        }

        if (!GA_is_continuous && SeqNo - TA_tracker == TR_length - GA_shift_front) {
            GA_tracker = SeqNo;

            // Calculate baseline average
            baseline_average /= (TR_length - GA_shift_back - TA_length - GA_shift_front);
            baseline_reset = true;
        }

        // Check if TA, TR, and GA are in progress. GA is an extension of TA that covers the ramp up and down of the gradient.
        TA_in_progress = TA_tracker > 0 && TA_tracker <= SeqNo && SeqNo < TA_tracker + TA_length;
        TR_in_progress = TA_tracker > 0 && TA_tracker <= SeqNo && SeqNo < TA_tracker + TR_length;
        GA_in_progress = GA_tracker > 0 && GA_tracker <= SeqNo && SeqNo < GA_tracker + GA_shift_front + TA_length + GA_shift_back;

        // Gradient artifact correction
        if (Apply_GACorr && GA_in_progress) {

            int temp_index = SeqNo - GA_tracker;

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

        // Apply baseline correction
        if (Apply_baseline && !GA_is_continuous) {

            if (GA_in_progress) {
                GA_sum += processing_sample_vector;
                processing_sample_vector += baseline_average - GA_average;
            } else if (TR_in_progress) {
                if (baseline_reset) {
                    GA_average = GA_sum / (GA_shift_front + TA_length + GA_shift_back);
                    GA_sum.setZero();
                    baseline_average.setZero();
                    baseline_reset = false;
                }
                baseline_average += processing_sample_vector;
            }
        }

        // Apply a geometric sum correction
        if (Apply_geometric_sum) {
            baseline_correction = (1 - baseline_update_rate) * baseline_correction + baseline_update_rate * processing_sample_vector;
            processing_sample_vector -= baseline_correction;
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
        if (getTriggerEnableStatus() && shouldTrigger(SeqNo) && checkTimeLimit() && !TA_in_progress) {
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
        
        // Samples, timestamp, triggers, and buffer index update
        sample_buffer_.col(current_data_index_) = processing_sample_vector;
        time_stamp_buffer_(current_data_index_) = time_stamp;
        trigger_buffer_A(current_data_index_) = trigger_A;
        trigger_buffer_B(current_data_index_) = trigger_B;
        
        // Save ROI means at the current index
        if (ROI_fMRI_means.size() > 0 && ROI_fMRI_means.size() == ROI_means_save.rows()) {
            ROI_means_save.col(current_data_index_) = ROI_fMRI_means;
        }
        
        current_sequence_number_ = SeqNo;
        current_data_index_ = (current_data_index_ + 1) % buffer_capacity_;

        if (current_data_index_ == 0) {
            sample_buffer_save.row(sample_buffer_save_index) = sample_buffer_.row(4);
            sample_buffer_save_index++;
            std::cout << "Row saved " << sample_buffer_save_index << std::endl;
        } else if (SeqNo >= 1500000 && data_saved == false) {
            writeMatrixdToCSV("sample_buffer_save.csv", sample_buffer_save);
            data_saved = true;
        }

    } catch (const std::exception& e) {
        std::cerr << "Datahandler exception: " << e.what() << '\n';
        std::cerr << boost::stacktrace::stacktrace();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;  // Explicitly specify duration type
    if(SeqNo > 100000) {
        total_addData_time += duration;
        min_addData_time = std::min(min_addData_time, duration);
        max_addData_time = std::max(max_addData_time, duration);
        addData_call_count++;
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

void dataHandler::updateSignalViewerData() {
    if (isReady()) {
        std::lock_guard<std::mutex> lock(this->dataMutex);
        if(sample_buffer_.rows() > 0) {
            int n_raw_channels = std::min(12, static_cast<int>(sample_buffer_.rows()));
            // Raw channels (0-11)
            for (int i = 0; i < n_raw_channels; i++) {
                emit channelDataUpdated(i, sample_buffer_.row(i), trigger_buffer_A, 
                                      trigger_buffer_B, trigger_buffer_out, time_stamp_buffer_, current_data_index_, channel_names_[i]);
            }

            // ROI channels (12+)
            for (int i = 0; i < ROI_means_save.rows(); i++) {
                emit channelDataUpdated(i + n_raw_channels, ROI_means_save.row(i), trigger_buffer_A, 
                                      trigger_buffer_B, trigger_buffer_out, time_stamp_buffer_, current_data_index_, ROI_names[i]);
            }
        }
    }
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

void dataHandler::savePreprocessingOutput(const Eigen::MatrixXd &output) {
    if (output.cols() == 0 || handler_state == WAITING_FOR_START) return;

    // Resize PreProcessing_output_save if needed
    if (PreProcessing_output_save.cols() != buffer_capacity_) {
        PreProcessing_output_save.resize(number_of_EEG_channels, buffer_capacity_);
        save_index_tracker = 0;
        current_data_index_ = 0;
    }

    int n_samples = output.cols();
    // Calculate the starting index to get samples from the middle
    int middle_start = (n_samples - n_samples) / 2;  // This will be 0 if n_samples == n_samples
    
    // Check if the new data would exceed buffer capacity
    if (save_index_tracker + n_samples > buffer_capacity_) {
        // Calculate how many samples can fit in the remaining space
        int remaining_space = buffer_capacity_ - save_index_tracker;
        
        // First part: Fill until the end of buffer, taking from middle of output
        PreProcessing_output_save.block(0, save_index_tracker, number_of_EEG_channels, remaining_space) = 
            output.block(0, middle_start, number_of_EEG_channels, remaining_space);
        
        // Second part: Wrap around and fill from the beginning
        int wrapped_samples = n_samples - remaining_space;
        if (wrapped_samples > 0) {
            PreProcessing_output_save.block(0, 0, number_of_EEG_channels, wrapped_samples) = 
                output.block(0, middle_start + remaining_space, number_of_EEG_channels, wrapped_samples);
        }
    } else {
        // Normal case: enough space to write all samples, taking from middle of output
        PreProcessing_output_save.block(0, save_index_tracker, number_of_EEG_channels, n_samples) = 
            output.block(0, middle_start, number_of_EEG_channels, n_samples);
    }
    
    // Update index with wrapping
    save_index_tracker = (save_index_tracker + n_samples) % buffer_capacity_;
}

void dataHandler::setROIMeans(const Eigen::VectorXd& means) {
    std::lock_guard<std::mutex> lock(this->dataMutex);
    if (ROI_means_save.rows() != means.size() || ROI_means_save.cols() != buffer_capacity_) {
        ROI_means_save.resize(means.size(), buffer_capacity_);
        ROI_means_save.setZero();
    }
    ROI_fMRI_means = means;
}

Eigen::VectorXd dataHandler::getROIMeanData(int roiIndex) {
    std::lock_guard<std::mutex> lock(this->dataMutex);
    if (roiIndex >= 0 && roiIndex < ROI_means_save.rows()) {
        return ROI_means_save.row(roiIndex);
    }
    return Eigen::VectorXd();  // Return empty vector if index is invalid
}
