#include "dataHandler.h"

// Resets the buffers for new channelcount and sampling rates
void dataHandler::reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate) {

    // Stop processing
    // processor_.stopProcessing();

    channel_count_ = channel_count;
    sampling_rate_ = sampling_rate;
    simulation_delivery_rate_ = simulation_delivery_rate;

    buffer_capacity_ = buffer_length_in_seconds_ * sampling_rate;
    current_data_index_ = 0;

    {
        std::lock_guard<std::mutex> lock(this->dataMutex);
        sample_buffer_ = Eigen::MatrixXd::Zero(channel_count, buffer_capacity_);
        time_stamp_buffer_ = Eigen::VectorXd::Zero(buffer_capacity_);
        trigger_buffer_ = Eigen::VectorXd::Zero(buffer_capacity_);
    }

    processing_sample_vector = Eigen::VectorXd::Zero(channel_count);

    GACorr_ = GACorrection(channel_count, GA_average_length, TA_length);

    baseline_average = Eigen::VectorXd(channel_count);

    handler_state = WAITING_FOR_STOP;

    RTfilter_.reset_filter(channel_count);
}
    
void dataHandler::reset_GACorr(int TA_length_input, int GA_average_length_input) {
    TA_length = TA_length_input; 
    GA_average_length = GA_average_length_input;
    GACorr_ = GACorrection(channel_count_, GA_average_length, TA_length);
    stimulation_tracker = 10000000;
}

int dataHandler::simulateData_sin() {
    auto startTime = std::chrono::high_resolution_clock::now();
    double time = 0.0;

    Eigen::VectorXd sample(channel_count_);
    const Eigen::VectorXd linspace_example = Eigen::VectorXd::LinSpaced(channel_count_, 0, channel_count_);

    int stimulation_interval = 10000;
    int stimulation_tracker = 0;
    while (true) { // Adjust this condition as needed
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = currentTime - startTime;
        double time_stamp = 0.0;
        int trigger = 0;
        
        // Check if it's time to generate the next sample
        if (elapsed.count() >= 1.0 / sampling_rate_) {
            startTime = currentTime;

            double SIN = 3.0 * M_PI * time;
            sample = linspace_example.unaryExpr([SIN](double x) { return std::sin(SIN * (5.0 + x)) * 20; });
            time_stamp = std::chrono::duration<double>(currentTime.time_since_epoch()).count();
            
            if (stimulation_tracker == stimulation_interval) {
                trigger = 1;
                stimulation_tracker = 0;
            } else {
                trigger = 0;
            }
            
            time += 1.0 / sampling_rate_;

            this->addData(sample, time_stamp, trigger, 0);

            stimulation_tracker++;
        }
    }

    return 0;
}

int dataHandler::simulateData_mat() {

    // Eigen::MatrixXd readMatFile("interl_eegfmri.mat");



    auto startTime = std::chrono::high_resolution_clock::now();
    double time = 0.0;

    Eigen::VectorXd sample(channel_count_);
    const Eigen::VectorXd linspace_example = Eigen::VectorXd::LinSpaced(channel_count_, 0, channel_count_);

    int stimulation_interval = 10000;
    int stimulation_tracker = 0;
    while (false) { // Adjust this condition as needed
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = currentTime - startTime;
        double time_stamp = 0.0;
        int trigger = 0;
        
        // Check if it's time to generate the next sample
        if (elapsed.count() >= 1.0 / sampling_rate_) {
            startTime = currentTime;

            double SIN = 5.0 * M_PI * time;
            sample = linspace_example.unaryExpr([SIN](double x) { return std::sin(SIN * (10.0 + x)); });
            time_stamp = std::chrono::duration<double>(currentTime.time_since_epoch()).count();
            
            if (stimulation_tracker >= stimulation_interval) {
                trigger = 1;
                stimulation_tracker = 0;
            } else {
                trigger = 0;
            }
            
            time += 1.0 / sampling_rate_;

            this->addData(sample, time_stamp, trigger, 0);

            stimulation_tracker++;
        }
    }

    return 0;
}

// For demonstration, print the size of a buffer (e.g., for channel 0)
void dataHandler::printBufferSize() {
    std::cout << "Buffer capacity: " << this->buffer_capacity_ << std::endl;
}

















// Buffer functions
// Add a signle sample to each channel
void dataHandler::addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger, const int &SeqNo) {

    if (trigger == 1) { 
        stimulation_tracker = 0;
        // GACorr_.printTemplate();
    }

    // Gradient artifact correction
    if (GACorr_running && stimulation_tracker < TA_length) {

        {
            std::lock_guard<std::mutex> lock(this->dataMutex);
            processing_sample_vector = samples - GACorr_.getTemplateCol(stimulation_tracker);
        }
        
        GACorr_.update_template(stimulation_tracker, samples);
        stimulation_tracker++;
    } else {

        std::lock_guard<std::mutex> lock(this->dataMutex);
        processing_sample_vector = samples;

    }

    // Baseline average
    if (Apply_baseline) {
        int baseline_end_index = (current_data_index_ - baseline_length + buffer_capacity_) % buffer_capacity_;
        baseline_average = baseline_average + (samples - sample_buffer_.col(baseline_end_index)) / baseline_length;
        processing_sample_vector = processing_sample_vector.cwiseQuotient(baseline_average);
    }

    // Filtering
    if (Apply_filter) {
        processing_sample_vector = RTfilter_.processSample(processing_sample_vector);
    }

    sample_buffer_.col(current_data_index_) = processing_sample_vector;

    // Triggering
    if (getTriggerEnableStatus() && shouldTrigger(SeqNo)) {
        send_trigger();
        removeTrigger(SeqNo);
    }

    // Timestamp, triggers and buffer index update
    time_stamp_buffer_(current_data_index_) = time_stamp;
    trigger_buffer_(current_data_index_) = trigger;
    current_sequence_number_ = SeqNo;
    current_data_index_ = (current_data_index_ + 1) % buffer_capacity_;
}

// Retrieves number_of_samples data points form all channels in chronological order
Eigen::MatrixXd dataHandler::returnLatestDataInOrder(int number_of_samples) {

    Eigen::MatrixXd output(channel_count_, number_of_samples);
    size_t channel_index = 0;

    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(number_of_samples, static_cast<int>(current_data_index_));
    int overflow = number_of_samples - fitToEnd;
    
    
    std::lock_guard<std::mutex> lock(this->dataMutex);
    if (fitToEnd > 0) {
        output.rightCols(fitToEnd) = sample_buffer_.middleCols(current_data_index_ - fitToEnd, fitToEnd);
    }
    
    if (overflow > 0) {
        output.leftCols(overflow) = sample_buffer_.middleCols(buffer_capacity_ - overflow, overflow);
    }

    return output;
}

// Retrieves number_of_samples data points form all channels in chronological order
int dataHandler::getLatestDataInOrder(Eigen::MatrixXd &output, int number_of_samples) {

    output.resize(channel_count_, number_of_samples);
    size_t channel_index = 0;

    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(number_of_samples, static_cast<int>(current_data_index_));
    int overflow = number_of_samples - fitToEnd;
    
    
    std::lock_guard<std::mutex> lock(this->dataMutex);
    if (fitToEnd > 0) {
        output.rightCols(fitToEnd) = sample_buffer_.middleCols(current_data_index_ - fitToEnd, fitToEnd);
    }
    
    if (overflow > 0) {
        output.leftCols(overflow) = sample_buffer_.middleCols(buffer_capacity_ - overflow, overflow);
    }

    return current_sequence_number_;
}

// Retrieves data from the specified channel in chronological order
Eigen::VectorXd dataHandler::getChannelDataInOrder(int channel_index, int downSamplingFactor) {

    if (downSamplingFactor < 1) throw std::invalid_argument("downSamplingFactor must be greater than 0");

    // Calculate the effective size after downsampling
    size_t downSampledSize = (buffer_capacity_ + downSamplingFactor - 1) / downSamplingFactor;
    Eigen::VectorXd downSampledData(downSampledSize);

    std::lock_guard<std::mutex> lock(this->dataMutex);
    for (size_t i = 0, j = 0; i < buffer_capacity_; i += downSamplingFactor, ++j) {
        // Calculate the index in the circular buffer accounting for wrap-around
        size_t index = (current_data_index_ + i) % buffer_capacity_;
        downSampledData(j) = sample_buffer_(channel_index, index);
    }

    return downSampledData;
}

// Retrieves data from the specified channels in chronological order
Eigen::MatrixXd dataHandler::getMultipleChannelDataInOrder(std::vector<int> channel_indices, int number_of_samples) {

    Eigen::MatrixXd outputData(channel_indices.size(), number_of_samples);
    size_t channel_index = 0;


    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(number_of_samples, (buffer_capacity_ - static_cast<int>(current_data_index_)));
    int overflow = number_of_samples - fitToEnd;
    
    
    std::lock_guard<std::mutex> lock(this->dataMutex);
    for (size_t j = 0; j < channel_indices.size(); ++j) {
        int channel_index = channel_indices[j];
        
        // Calculate the number of samples that fit to the end and the overflow
        int overflow = number_of_samples - fitToEnd;

        if (fitToEnd > 0) {
            outputData.row(j).head(fitToEnd) = sample_buffer_.block(channel_index, current_data_index_, 1, fitToEnd);
        }
        
        if (overflow > 0) {
            outputData.row(j).tail(overflow) = sample_buffer_.block(channel_index, 0, 1, overflow);
        }
    }

    return outputData;
}

// More efficient alternative to getMultipleChannelDataInOrder that can be used if the wanted channels are adjacent. Retrieves adjacent channels as a block
Eigen::MatrixXd dataHandler::getBlockChannelDataInOrder(int first_channel_index, int number_of_channels, int number_of_samples) {

    Eigen::MatrixXd outputData(number_of_channels, number_of_samples);
    size_t channel_index = 0;

    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(number_of_samples, static_cast<int>(current_data_index_));
    int overflow = number_of_samples - fitToEnd;
    
    
    std::lock_guard<std::mutex> lock(this->dataMutex);
    if (fitToEnd > 0) {
        outputData.rightCols(fitToEnd) = sample_buffer_.block(first_channel_index, current_data_index_ - fitToEnd, number_of_channels, fitToEnd);
    }
    
    if (overflow > 0) {
        outputData.leftCols(overflow) = sample_buffer_.block(first_channel_index, buffer_capacity_ - overflow, number_of_channels, overflow);
    }

    return outputData;
}

// Retrieves data from the time stamp channel in chronological order
Eigen::VectorXd dataHandler::getTimeStampsInOrder(int downSamplingFactor) {

    std::lock_guard<std::mutex> lock(this->dataMutex);

    if (downSamplingFactor == 0) throw std::invalid_argument("downSamplingFactor must be greater than 0");

    // Calculate the effective size after downsampling
    size_t downSampledSize = (buffer_capacity_ + downSamplingFactor - 1) / downSamplingFactor;
    Eigen::VectorXd downSampledData(downSampledSize);

    for (size_t i = 0, j = 0; i < buffer_capacity_; i += downSamplingFactor, ++j) {
        // Calculate the index in the circular buffer accounting for wrap-around
        size_t index = (current_data_index_ + i) % buffer_capacity_;
        downSampledData(j) = time_stamp_buffer_(index);
    }

    return downSampledData;
}

// Retrieves data from the trigger channel in chronological order
Eigen::VectorXd dataHandler::getTriggersInOrder(int downSamplingFactor) {

    std::lock_guard<std::mutex> lock(this->dataMutex);

    if (downSamplingFactor == 0) throw std::invalid_argument("downSamplingFactor must be greater than 0");

    // Calculate the effective size after downsampling
    size_t downSampledSize = (buffer_capacity_ + downSamplingFactor - 1) / downSamplingFactor;
    Eigen::VectorXd downSampledData(downSampledSize);

    for (size_t i = 0, j = 0; i < buffer_capacity_; i += downSamplingFactor, ++j) {
        // Calculate the index in the circular buffer accounting for wrap-around
        size_t index = (current_data_index_ + i) % buffer_capacity_;
        downSampledData(j) = trigger_buffer_(index);
    }

    return downSampledData;
}








// Trigger sending

int dataHandler::connectTriggerPort() {
    return magPro_.connectTriggerPort();
}

void dataHandler::send_trigger() {
    magPro_.trig();
}

void dataHandler::set_enable(bool status) {
    magPro_.set_enable(status);
    triggerEnableState = status;
}

void dataHandler::set_amplitude(int amplitude) {
    magPro_.set_amplitude(amplitude);
}
