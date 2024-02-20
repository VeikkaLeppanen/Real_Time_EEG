#include "dataHandler.h"

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
        if (elapsed.count() >= 1.0 / delivery_rate_) {
            startTime = currentTime;

            // Single sample per packet
            if (sample_packet_size_ == 1) {
                double SIN = 2.0 * M_PI * time;
                sample = linspace_example.unaryExpr([SIN](double x) { return std::sin(SIN * (1.0 + x)); });
                sample(sample.size() - 1) = std::chrono::duration<double>(currentTime.time_since_epoch()).count();
                time += 1.0 / sampling_rate_;

                {
                std::unique_lock<std::mutex> mlock(dataMutex, std::try_to_lock);                
                if (mlock) { this->addData(sample); }
                }

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

                {
                std::unique_lock<std::mutex> mlock(dataMutex, std::try_to_lock);                
                if (mlock) { this->addData(sample_packet); }
                }
            }
        }
    }
    return 0;
}

// Add samples to channels. Last row in the shortBuffer_ is for the time_stamp
template<typename Derived>
void dataHandler::addData(const Eigen::MatrixBase<Derived> &sample_packet) {
    if (sample_packet.rows() != channel_count_ + 1) {
        std::cerr << "Invalid sample dimensions: " << sample_packet.rows() << " expected: " << channel_count_ + 1 << std::endl;
        return;
    }
    shortBuffer_.addSamples(sample_packet.derived());
}

// For demonstration, print the size of a buffer (e.g., for channel 0)
void dataHandler::printBufferSize(int channel) {
    if (channel < 0 || channel > channel_count_) {
        std::cerr << "Invalid channel index" << std::endl;
        return;
    }
    // Assuming circularEigenBuffer has a method getCapacity()
    std::cout << "Buffer capacity: " << shortBuffer_.getCapacity() << std::endl;
}

