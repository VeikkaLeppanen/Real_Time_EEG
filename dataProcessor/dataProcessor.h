#ifndef DATAPROCESSOR_H
#define DATAPROCESSOR_H

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <Eigen/Dense>

class dataProcessor {
public:
    dataProcessor(size_t threadCount)
        : running(true), current_data_index_(0), last_processed_index_(0) {
        workers.reserve(threadCount);
        for (size_t i = 0; i < threadCount; ++i) {
            workers.emplace_back(&dataProcessor::processData, this);
        }
    }

    ~dataProcessor() {
        // Stop all processing
        stopProcessing();
        for (auto& worker : workers) {
            if (worker.joinable())
                worker.join();
        }
    }

    void newData(const Eigen::VectorXd &samples) {
        std::unique_lock<std::mutex> lock(data_mutex);
        sample_buffer_.col(current_data_index_) = samples;
        current_data_index_ = (current_data_index_ + 1) % buffer_capacity_;
        lock.unlock();
        cv.notify_one(); 
    }

    void reset_processor(int channel_count, int buffer_capacity) {
        std::lock_guard<std::mutex> lock(data_mutex);
        channel_count_ = channel_count;
        buffer_capacity_ = buffer_capacity;
        sample_buffer_ = Eigen::MatrixXd::Zero(channel_count, buffer_capacity);
        current_data_index_ = 0;
        last_processed_index_ = 0;
    }

    void startProcessing() {
        std::unique_lock<std::mutex> lock(data_mutex);
        processing = true;
        cv.notify_all();
    }

    void stopProcessing() {
        std::unique_lock<std::mutex> lock(data_mutex);
        processing = false;
    }

    void stopAll() {
        std::unique_lock<std::mutex> lock(data_mutex);
        running = false;
        cv.notify_all();
    }

private:
    std::vector<std::thread> workers;
    std::mutex data_mutex;
    std::condition_variable cv;
    bool running;
    bool processing;

    Eigen::MatrixXd sample_buffer_;
    int buffer_capacity_;
    int channel_count_;
    size_t current_data_index_;
    size_t last_processed_index_;

    Eigen::MatrixXd getBlockChannelDataInOrder(int first_channel_index, int number_of_channels, int number_of_samples);

    Eigen::MatrixXd downsample(const Eigen::MatrixXd& input, int factor);
    Eigen::MatrixXd delayEmbed(const Eigen::MatrixXd& CWL, int delay);
    Eigen::MatrixXd removeBCG(const Eigen::MatrixXd& EEG, Eigen::MatrixXd CWL, int delay);

    void processData() {
        do {
            while (processing && running) {
                std::unique_lock<std::mutex> lock(data_mutex);
                std::cout << "Thread waiting...\n";  // Debug statement
                if (current_data_index_ == last_processed_index_) {
                    cv.wait(lock, [this]{ return current_data_index_ != last_processed_index_ || !processing || !running; });
                    if (!processing || !running) break;
                }

                // Suppose your actual processing requires current data
                if (current_data_index_ == last_processed_index_) {
                    std::cout << "No new data to process.\n";  // Debug statement
                    continue;  // Skip processing if no new data
                }

                // Copy data from wanted EEG channels and CWL channels
                // Make sure these functions are thread-safe and do not block indefinitely
                Eigen::MatrixXd EEG = getBlockChannelDataInOrder(0, 4, 10000);
                Eigen::MatrixXd CWL = getBlockChannelDataInOrder(8, 4, 10000);

                last_processed_index_ = current_data_index_;  // Update processed index
                lock.unlock();  // Unlock early to allow other threads to access the buffer

                // Process the data
                process(EEG, CWL);
            }
        } while (running);
    }

    void process(Eigen::MatrixXd EEG, Eigen::MatrixXd CWL);
};


#endif // DATAPROCESSOR_H
