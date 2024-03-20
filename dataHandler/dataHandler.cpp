#include "dataHandler.h"

// Resets the buffers for new channelcount and sampling rates
void dataHandler::reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate) {

    this->channel_count_ = channel_count;
    this->sampling_rate_ = sampling_rate;
    this->simulation_delivery_rate_ = simulation_delivery_rate;

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

            this->addData(sample, time_stamp, trigger);

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

            this->addData(sample, time_stamp, trigger);

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
void dataHandler::addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger) {

    if (trigger == 1) { 
        stimulation_tracker = 0;
        // GACorr_.printTemplate();
    }

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

    time_stamp_buffer_(current_data_index_) = time_stamp;
    trigger_buffer_(current_data_index_) = trigger;
    current_data_index_ = (current_data_index_ + 1) % buffer_capacity_;
}

// Retrieves data form all channels in chronological order
Eigen::MatrixXd dataHandler::getDataInOrder(int downSamplingFactor) {
    
    // LOWPASS filter 120Hz   BANDPASS 0.33-125Hz
    // DONWSAMPLINGFACTOR 10

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

// Retrieves data from the time stamp channel in chronological order
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

// Retrieves data from the trigger channel in chronological order
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






// Converts .mat files into an eigen matrix
// Eigen::MatrixXd dataHandler::readMatFile(const std::string& fileName) {
//     mat_t *matfp = Mat_Open(fileName.c_str(), MAT_ACC_RDONLY);
//     if (matfp == nullptr) {
//         throw std::runtime_error("Error opening MAT file.");
//     }

//     // Replace "variableName" with the name of your variable in the MAT file
//     matvar_t *matvar = Mat_VarRead(matfp, "interleavedEegfmri");
//     if (matvar == nullptr) {
//         Mat_Close(matfp);
//         throw std::runtime_error("Error reading variable from MAT file.");
//     }

//     if (matvar->rank != 2 || matvar->data_type != MAT_T_DOUBLE) {
//         Mat_VarFree(matvar);
//         Mat_Close(matfp);
//         throw std::runtime_error("Variable must be a 2D double array.");
//     }

//     size_t rows = matvar->dims[0];
//     size_t cols = matvar->dims[1];
//     double* data = static_cast<double*>(matvar->data);

//     // Transfer data to Eigen
//     Eigen::MatrixXd matrix(rows, cols);
//     for (size_t i = 0; i < rows; ++i) {
//         for (size_t j = 0; j < cols; ++j) {
//             matrix(i, j) = data[i + j * rows]; // Column-major order in MATLAB
//         }
//     }

//     Mat_VarFree(matvar);
//     Mat_Close(matfp);

//     return matrix;
// }
