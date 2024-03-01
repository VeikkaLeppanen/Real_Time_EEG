#include "dataHandler.h"

// Resets the buffers for new channelcount and sampling rates
void dataHandler::reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate) {

    this->channel_count_ = channel_count;
    this->sampling_rate_ = sampling_rate;
    this->simulation_delivery_rate_ = simulation_delivery_rate;
    this->sample_packet_size_ = sampling_rate_ / simulation_delivery_rate_;

    this->short_buffer_capacity_ = short_buffer_length_in_seconds_ * sampling_rate;
    this->long_buffer_capacity_ = long_buffer_length_in_seconds_ * sampling_rate;

    {
        std::lock_guard<std::mutex> (this->dataMutex);
        this->short_buffer_ = circularEigenBuffer(channel_count + 1, short_buffer_capacity_);
        this->long_buffer_ = circularEigenBuffer(channel_count + 1, long_buffer_capacity_);
    }
}

int dataHandler::simulateData() {
    auto startTime = std::chrono::high_resolution_clock::now();
    double time = 0.0;

    Eigen::MatrixXd sample_packet(channel_count_ + 1, sample_packet_size_);
    Eigen::VectorXd sample(channel_count_ + 1);
    const Eigen::VectorXd linspace_example = Eigen::VectorXd::LinSpaced(channel_count_ + 1, 0, channel_count_);

    while (true) { // Adjust this condition as needed
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = currentTime - startTime;
        
        // Check if it's time to generate the next sample
        if (elapsed.count() >= 1.0 / simulation_delivery_rate_) {
            startTime = currentTime;

            // Single sample per packet
            if (sample_packet_size_ == 1) {
                double SIN = 2.0 * M_PI * time;
                sample = linspace_example.unaryExpr([SIN](double x) { return std::sin(SIN * (1.0 + x)); });
                sample(sample.size() - 1) = std::chrono::duration<double>(currentTime.time_since_epoch()).count();
                time += 1.0 / sampling_rate_;

                this->addData(sample);

            // Multiple samples per packet
            } else {
                for (int i = 0; i < sample_packet_size_; i++) {
                    double SIN = 2.0 * M_PI * time;
                    sample = linspace_example.unaryExpr([SIN](double x) { return std::sin(SIN * (1.0 + x)); });
                    sample_packet.col(i) = sample;
                    time += 1.0 / sampling_rate_;
                }

                double time_stamp = std::chrono::duration<double>(currentTime.time_since_epoch()).count();
                sample_packet.row(sample.size() - 1).fill(time_stamp);

                this->addData(sample_packet);
            }
        }
    }

    return 0;
}

// Add samples to channels. Last row in the short_buffer_ is for the time_stamp
template<typename Derived>
void dataHandler::addData(const Eigen::MatrixBase<Derived> &sample_packet) {
    if (sample_packet.rows() != channel_count_ + 1) {
        std::cerr << "Invalid sample dimensions: " << sample_packet.rows() << " expected: " << channel_count_ + 1 << std::endl;
        return;
    }

    {
        std::lock_guard<std::mutex> guard(this->dataMutex);

        // Save samples to the short buffer
        size_t old_index = short_buffer_.getCurrentIndex();
        short_buffer_.addSamples(sample_packet.derived());

        // If the short buffer is filled copy it to the long buffer
        if (useLongBuffer && short_buffer_.getCurrentIndex() < old_index) {
            long_buffer_.addSamples(short_buffer_.getDataInOrder());
        }
    }
}


// Retrieves data from the specified channel in chronological order
Eigen::VectorXd dataHandler::getChannelDataInOrder(int channel_index, int downSamplingFactor) {
    {
        std::lock_guard<std::mutex> guard(this->dataMutex);
        return this->short_buffer_.getChannelDataInOrder(channel_index, downSamplingFactor);
    }
}

// Retrieves data form all channels in chronological order
Eigen::MatrixXd dataHandler::getDataInOrder(int downSamplingFactor) {
    {
        std::lock_guard<std::mutex> guard(this->dataMutex);
        return this->short_buffer_.getDataInOrder(downSamplingFactor);
    }
}

// For demonstration, print the size of a buffer (e.g., for channel 0)
void dataHandler::printBufferSize() {
    std::cout << "Buffer capacity: " << short_buffer_.getCapacity() << std::endl;
}

