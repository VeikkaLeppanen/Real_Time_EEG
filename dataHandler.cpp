#include "dataHandler.h"

int dataHandler::simulateData(int sampling_rate) {
    auto startTime = std::chrono::high_resolution_clock::now();
    double time = 0.0;

    while (true) { // Adjust this condition as needed
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = currentTime - startTime;

        // Check if it's time to generate the next sample
        if (elapsed.count() >= 1.0 / sampling_rate) {
            startTime = currentTime;

            Eigen::VectorXd samples = Eigen::VectorXd::LinSpaced(channel_count_ + 1, 0, channel_count_);
            double SIN = 2.0 * M_PI * time;
            samples = samples.unaryExpr([SIN](double x) { return std::sin(SIN * (1.0 + x)); });

            {
            std::unique_lock<std::mutex> mlock(dataMutex, std::try_to_lock);                
            if (mlock) { this->addData(samples, std::chrono::duration<double>(currentTime.time_since_epoch()).count()); }
            }

            time += 1.0 / sampling_rate;
        }
    }
    return 0;
}

void dataHandler::addData(Eigen::VectorXd samples, double time_stamp) {
    if (samples.size() != channel_count_ + 1) {
        std::cerr << "Invalid sample dimensions: " << samples.size() << " expected: " << channel_count_ + 1 << std::endl;
        return;
    }
    samples(samples.size() - 1) = time_stamp;
    dataBuffer.addSamples(samples);
}


// TODO: implement correctly
void dataHandler::addData(Eigen::MatrixXd samples, double time_stamp) {
    if (samples.rows() != channel_count_) {
        std::cerr << "Invalid sample dimensions: " << samples.cols() << " expected: " << channel_count_ << std::endl;
        return;
    }
    dataBuffer.addSamples(samples);
}

// For demonstration, print the size of a buffer (e.g., for channel 0)
void dataHandler::printBufferSize(int channel) {
    if (channel < 0 || channel > channel_count_) {
        std::cerr << "Invalid channel index" << std::endl;
        return;
    }
    // Assuming circularEigenBuffer has a method getCapacity()
    std::cout << "Buffer capacity: " << dataBuffer.getCapacity() << std::endl;
}

