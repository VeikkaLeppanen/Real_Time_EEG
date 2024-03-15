#include "dataHandler.h"

// Resets the buffers for new channelcount and sampling rates
void dataHandler::reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate) {

    this->channel_count_ = channel_count;
    this->sampling_rate_ = sampling_rate;
    this->simulation_delivery_rate_ = simulation_delivery_rate;
    this->sample_packet_size_ = sampling_rate_ / simulation_delivery_rate_;

    this->buffer_capacity_ = buffer_length_in_seconds_ * sampling_rate;
    this->current_data_index_ = 0;

    {
        std::lock_guard<std::mutex> (this->dataMutex);
        this->sample_buffer_ = Eigen::MatrixXd::Zero(channel_count, buffer_capacity_);
        this->time_stamp_buffer_ = Eigen::VectorXd::Zero(buffer_capacity_);
        this->trigger_buffer_ = Eigen::VectorXd::Zero(buffer_capacity_);
    }

    this->GACorr_ = GACorrection(channel_count, this->GA_average_length, this->TA_length);
}

int dataHandler::simulateData() {
    auto startTime = std::chrono::high_resolution_clock::now();
    double time = 0.0;

    Eigen::MatrixXd sample_packet(channel_count_, sample_packet_size_);
    Eigen::VectorXd sample(channel_count_);
    const Eigen::VectorXd linspace_example = Eigen::VectorXd::LinSpaced(channel_count_, 0, channel_count_);

    int stimulation_interval = 4000;
    int stimulation_tracker = 0;
    while (true) { // Adjust this condition as needed
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = currentTime - startTime;
        Eigen::VectorXd time_stamps = Eigen::VectorXd::Zero(sample_packet_size_);
        Eigen::VectorXd triggers = Eigen::VectorXd::Zero(sample_packet_size_);
        
        // Check if it's time to generate the next sample
        if (elapsed.count() >= 1.0 / simulation_delivery_rate_) {
            startTime = currentTime;

            // Single sample per packet
            if (sample_packet_size_ == 1) {
                double SIN = 5.0 * M_PI * time;
                sample = linspace_example.unaryExpr([SIN](double x) { return std::sin(SIN * (10.0 + x)); });
                time_stamps.fill(std::chrono::duration<double>(currentTime.time_since_epoch()).count());
                
                if (stimulation_tracker == stimulation_interval) {
                    triggers.fill(1);
                    stimulation_tracker = 0;
                } else {
                    triggers.fill(0);
                }
                
                time += 1.0 / sampling_rate_;

                this->addData(sample, time_stamps, triggers);

            // Multiple samples per packet
            } else {
                for (int i = 0; i < sample_packet_size_; i++) {
                    double SIN = 2.0 * M_PI * time;
                    sample = linspace_example.unaryExpr([SIN](double x) { return std::sin(SIN * (1.0 + x)); });
                    sample_packet.col(i) = sample;
                    time += 1.0 / sampling_rate_;
                }

                time_stamps.fill(std::chrono::duration<double>(currentTime.time_since_epoch()).count());
                triggers.fill(0);

                this->addData(sample_packet, time_stamps, triggers);
            }
            stimulation_tracker++;
        }
    }

    return 0;
}

// TODO: This needs to be implemented to accept sample packets with multiple bundles.
// Add samples to channels. First rows of sample packet are the data and the last two are stimulation flag and time stamp.
// void dataHandler::addDataGACorr(const Eigen::VectorXd &sample_packet) {
//     if (sample_packet.rows() != channel_count_) {
//         std::cerr << "Invalid sample dimensions: " << sample_packet.rows() << " expected: " << channel_count_<< std::endl;
//         return;
//     }

//     // Check for stimulation flag
//     if (sample_packet(sample_packet.size() - 2) == 1) { this->stimulation_tracker = 0; }

//     if (this->stimulation_tracker < this->TA_length) {
        
//         Eigen::VectorXd samples_corrected = sample_packet - GACorr_.getTemplateCol(this->stimulation_tracker);
//         std::lock_guard<std::mutex> guard(this->dataMutex);
//         sample_buffer_.addSamples(samples_corrected);

//         this->stimulation_tracker++;
//     } else { // Add data to buffer
//         std::lock_guard<std::mutex> guard(this->dataMutex);
//         sample_buffer_.addSamples(sample_packet.derived());
//     }
// }

// template<typename Derived>
// void dataHandler::addData(const Eigen::MatrixBase<Derived> &sample_packet, const Eigen::VectorXd &time_stamps, const Eigen::VectorXd &triggers) {
//     if (sample_packet.rows() != channel_count_) {
//         std::cerr << "Invalid sample dimensions: " << sample_packet.rows() << " expected: " << channel_count_ << std::endl;
//         return;
//     }

//     { // Add data to buffer
//         std::lock_guard<std::mutex> guard(this->dataMutex);
//         sample_buffer_.addSamples(sample_packet.derived());
//         time_stamp_buffer_.addSamples(time_stamp.transpose());
//         trigger_buffer_.addSamples(triggers.transpose());
//     }
// }

// Retrieves data from the specified channel in chronological order
// Eigen::VectorXd dataHandler::getChannelDataInOrder(int channel_index, int downSamplingFactor) {
//     {
//         std::lock_guard<std::mutex> guard(this->dataMutex);
//         return this->sample_buffer_.getChannelDataInOrder(channel_index, downSamplingFactor);
//     }
// }

// // Retrieves data form all channels in chronological order
// Eigen::MatrixXd dataHandler::getDataInOrder(int downSamplingFactor) {
//     {
//         std::lock_guard<std::mutex> guard(this->dataMutex);
//         return this->sample_buffer_.getDataInOrder(downSamplingFactor);
//     }
// }

// For demonstration, print the size of a buffer (e.g., for channel 0)
void dataHandler::printBufferSize() {
    std::cout << "Buffer capacity: " << this->buffer_capacity_ << std::endl;
}

















// Buffer functions
// Add a signle sample to each channel
void dataHandler::addData(const Eigen::VectorXd &samples, const Eigen::VectorXd &time_stamps, const Eigen::VectorXd &triggers) {

    if (triggers(0) == 1.0) { stimulation_tracker = 0; }

    if (stimulation_tracker < TA_length) {

        {
            std::lock_guard<std::mutex> (this->dataMutex);
            sample_buffer_.col(current_data_index_) = samples - GACorr_.getTemplateCol(stimulation_tracker);
        }
        
        GACorr_.update_template(stimulation_tracker, samples);
        stimulation_tracker++;
    } else {

        std::lock_guard<std::mutex> (this->dataMutex);
        sample_buffer_.col(current_data_index_) = samples;

    }

    time_stamp_buffer_(current_data_index_) = time_stamps(0);
    trigger_buffer_(current_data_index_) = triggers(0);
    current_data_index_ = (current_data_index_ + 1) % buffer_capacity_;
}

// Add multiple samples to each channel
void dataHandler::addData(const Eigen::MatrixXd &samples, const Eigen::VectorXd &time_stamps, const Eigen::VectorXd &triggers) {

    std::lock_guard<std::mutex> (this->dataMutex);

    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(samples.cols(), static_cast<Eigen::Index>(buffer_capacity_ - current_data_index_));
    
    // Assign samples that fit to the end
    if (fitToEnd > 0) {
        sample_buffer_.middleCols(current_data_index_, fitToEnd) = samples.leftCols(fitToEnd);
    }
    
    // Calculate remaining samples that need to wrap around
    int overflow = samples.cols() - fitToEnd;
    
    // If there are any samples that didn't fit, wrap them around to the beginning
    if (overflow > 0) {
        sample_buffer_.leftCols(overflow) = samples.rightCols(overflow);
    }

    current_data_index_ = (current_data_index_ + samples.cols()) % buffer_capacity_;
}

// Retrieve the sample at a specific index in the rolling buffer
// Index 0 refers to the oldest sample, and capacity-1 refers to the newest sample.
// double dataHandler::getSample(int channel, int index) const {
//     if (index >= row_capacity_) throw std::out_of_range("Buffer index out of bounds");

//     // Calculate the actual index based on currentIndex
//     int actualIndex = (currentIndex_ + index) % row_capacity_;
//     return data_(channel, actualIndex);
// }

// Retrieves data form all channels in chronological order
Eigen::MatrixXd dataHandler::getDataInOrder(int downSamplingFactor) {

    std::lock_guard<std::mutex> (this->dataMutex);

    if (downSamplingFactor == 0) throw std::invalid_argument("downSamplingFactor must be greater than 0");

    // Calculate the effective size after downsampling
    size_t downSampledSize = (buffer_capacity_ + downSamplingFactor - 1) / downSamplingFactor;
    Eigen::MatrixXd downSampledData(sample_buffer_.rows(), downSampledSize);

    for (size_t i = 0, j = 0; i < buffer_capacity_; i += downSamplingFactor, ++j) {
        // Calculate the index in the circular buffer accounting for wrap-around
        size_t index = (current_data_index_ + i) % buffer_capacity_;
        downSampledData.col(j) = sample_buffer_.col(index);
    }

    return downSampledData;
}

// Retrieves data from the specified channel in chronological order
Eigen::VectorXd dataHandler::getChannelDataInOrder(int channel_index, int downSamplingFactor) {

    std::lock_guard<std::mutex> (this->dataMutex);

    if (downSamplingFactor == 0) throw std::invalid_argument("downSamplingFactor must be greater than 0");

    // Calculate the effective size after downsampling
    size_t downSampledSize = (buffer_capacity_ + downSamplingFactor - 1) / downSamplingFactor;
    Eigen::VectorXd downSampledData(downSampledSize);

    for (size_t i = 0, j = 0; i < buffer_capacity_; i += downSamplingFactor, ++j) {
        // Calculate the index in the circular buffer accounting for wrap-around
        size_t index = (current_data_index_ + i) % buffer_capacity_;
        downSampledData(j) = sample_buffer_(channel_index, index);
    }

    return downSampledData;
}

// Retrieves data from the specified channel in chronological order
Eigen::VectorXd dataHandler::getTimeStampsInOrder(int downSamplingFactor) {

    std::lock_guard<std::mutex> (this->dataMutex);

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

// Retrieves data from the specified channel in chronological order
Eigen::VectorXd dataHandler::getTriggersInOrder(int downSamplingFactor) {

    std::lock_guard<std::mutex> (this->dataMutex);

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