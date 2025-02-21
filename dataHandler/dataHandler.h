#ifndef DATAHANDLER_H
#define DATAHANDLER_H

#include <mutex>
#include <condition_variable>
#include <queue>

#include <cstddef> // For size_t
#include <iostream>
#include <chrono>
#include <array>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>
#include <LabJackM.h>
#include <QObject>
#include <QTimer>
#include <QDebug>

#include "EEG/preprocessing/GACorrection.h"
#include "devices/TMS/magPro/magPro.h"
#include "../EEG/preprocessing/preprocessingFunctions.h"
#include "../utils/utilityFunctions.h"
#include "devices/EEG/eeg_bridge/eeg_bridge.h"
#include <boost/stacktrace.hpp>

// In case the bind fails, use the following commands to find the PID and kill the process:
// lsof -i :8080
// kill -9 PID

enum HandlerState {
  WAITING_FOR_START,
  WAITING_FOR_STOP
};

enum TMSConnectionType {
    COM,
    TTL
};

class dataHandler : public QObject {
    Q_OBJECT
public:
    // Default constructor
    explicit dataHandler(QObject *parent = nullptr)  // dataProcessor &processor
                :   QObject(parent),
                    channel_count_(0),
                    sampling_rate_(0),
                    simulation_delivery_rate_(0),
                    GACorr_(GACorrection(0, 0, 0)),
                    magPro_3G()
    {
        // Move timer creation and setup to the main thread using moveToThread
        QMetaObject::invokeMethod(this, [this]() {
            qDebug() << "Creating dataHandler timer";
            timer = new QTimer(this);
            connect(timer, &QTimer::timeout, this, &dataHandler::updateSignalViewerData, Qt::QueuedConnection);
            timer->start(16);
        }, Qt::QueuedConnection);
    }

    // Reset functions
    void reset_handler(int channel_count, int sampling_rate, int simulation_delivery_rate);
    void reset_handler(int channel_count, int sampling_rate) { reset_handler(channel_count, sampling_rate, sampling_rate); }

    void setHandlerState(HandlerState state) { handler_state = state; }
    bool isReady() { return (handler_state == WAITING_FOR_STOP); }

    // Data handling
    void addData(const Eigen::VectorXd &samples, const double &time_stamp, const int &trigger_A, const int &trigger_B, const int &SeqNo);

    int getLatestSequenceNumber() { return current_sequence_number_; }
    int getLatestDataAndTriggers(Eigen::MatrixXd &output, 
                                 Eigen::VectorXi &triggers_A, 
                                 Eigen::VectorXi &triggers_B, 
                                 Eigen::VectorXi &triggers_out, 
                                 Eigen::VectorXd &time_stamps, 
                                             int number_of_samples);

    int get_buffer_capacity() { return buffer_capacity_; }
    int get_buffer_length_in_seconds() { return buffer_length_in_seconds_; }
    int get_channel_count() { return channel_count_; }
    int get_ROI_channel_count() { return ROI_fMRI_means.size(); }
    int get_current_data_index() { return current_data_index_; }
    int getSamplingRate() { return sampling_rate_; }

    int getPreprocessingOutput(Eigen::MatrixXd &output, 
                                 Eigen::VectorXi &triggers_A, 
                                 Eigen::VectorXi &triggers_B, 
                                 Eigen::VectorXi &triggers_out, 
                                 Eigen::VectorXd &time_stamps, 
                                             int &number_of_samples) 
    { 
        std::lock_guard<std::mutex> lock(this->dataMutex);
        output = preprocessing_output;
        triggers_A = preprocessing_triggers_A;
        triggers_B = preprocessing_triggers_B;
        triggers_out = preprocessing_triggers_out;
        time_stamps = preprocessing_time_stamps;
        number_of_samples = preprocessing_number_of_samples;
        return processing_sequence_number;
    }

    void setTriggerSource(uint16_t source) { triggerSource = source; }
    void setSourceChannels(std::vector<uint16_t> SourceChannels) {
        source_channels_.resize(SourceChannels.size());
        for (size_t i = 0; i < SourceChannels.size(); ++i) {
            source_channels_(static_cast<int>(i)) = static_cast<int>(SourceChannels[i]);
        }
    }
    Eigen::VectorXi getSourceChannels() const { return source_channels_; }
    
    void setChannelNames(std::vector<std::string> channel_names) { 
        channel_names_ = channel_names;
        channel_names_set = true;
        emit channelNamesUpdated(channel_names_);
    }
    bool channelNamesSet() { return channel_names_set; }
    std::vector<std::string> getChannelNames() { return channel_names_; }

    // Gradient artifact correction
    void GACorr_off() { Apply_GACorr = false; }
    void GACorr_on() {
        int TA_tracker = 10000000;
        Apply_GACorr = true; 
    }

    void reset_GACorr(int TA_length_input, int GA_average_length_input);
    void reset_GACorr_tracker() { 
        TA_tracker = 10000000;
        GACorr_.reset_index();
    }

    void set_TR_length(int TR_length) { TR_length = TR_length; }

    int get_TA_length() { return TA_length; }
    int get_TR_length() { return TR_length; }
    int get_GA_average_length() { return GA_average_length; }
    bool getGAState() { return Apply_GACorr; }
    
    // Filtering
    void setFilterState(bool state) { Apply_filter = state; }
    bool getFilterState() { return Apply_filter; }

    // Baseline
    void setBaselineState(bool state) { Apply_baseline = state; }
    bool getBaselineState() { return Apply_baseline; }

    // Geometric sum correction
    void setGeometricSumState(bool state) { Apply_geometric_sum = state; }
    bool getGeometricSumState() { return Apply_geometric_sum; }

    // Triggering
    void setTriggerConnectStatus(bool value) { triggerPortState = value; }
    bool getTriggerConnectStatus() { return triggerPortState; }

    void setTriggerEnableStatus(bool value) { triggerEnableState = value; }
    bool getTriggerEnableStatus() { return triggerEnableState; }

    void setTriggerTimeLimit(int value) { time_limit = std::max(min_time_limit, std::min(max_time_limit, value)); }
    double getTriggerTimeLimit() { return time_limit; }
    bool checkTimeLimit() {
        // Calculate the time difference and compare it to the time_limit
        auto duration_since_trigger = std::chrono::system_clock::now() - latest_trigger_time;
        auto duration_limit = std::chrono::milliseconds(time_limit);

        bool enough_time_passed = duration_since_trigger > duration_limit;
        return enough_time_passed;
    }

    void insertTrigger(int seqNum) {
        std::lock_guard<std::mutex> lock(triggerMutex);
        triggerSet.insert(seqNum);
    }

    bool shouldTrigger(int seqNum) {
        std::lock_guard<std::mutex> lock(triggerMutex);
        return triggerSet.find(seqNum) != triggerSet.end();
    }

    void removeTrigger(int seqNum) {
        std::lock_guard<std::mutex> lock(triggerMutex);
        triggerSet.erase(seqNum);
    }

    void setTMSConnectionType(TMSConnectionType type) { TMS_connectionType = type; }

    // MAGPRO FUNCTIONS
    int connectTriggerPort();
    void send_trigger();
    void set_enable(bool status);
    void set_amplitude(int amplitude);
    void magPro_set_mode(int mode = 0, int direction = 0, int waveform = 1, int burst_pulses = 5, float ipi = 1, float ba_ratio = 1.0, bool delay = true);
    void magPro_request_mode_info();
    void get_mode_info(int &mode, int &direction, int &waveform, int &burst_pulses, float &ipi, float &ba_ratio, bool &enabled);
    void save_seqnum_list() { 
        writeMatrixiToCSV("trigger_seqNum_list.csv", vectorToColumnMatrixi(seqNum_list)); 
        qDebug() << "Trigger list size:" << seqNum_list.size();

        // Print timing statistics for addData
        if (addData_call_count > 0) {
            double avg_time = total_addData_time.count() / addData_call_count;
            qDebug() << "addData statistics:"
                     << "\nTotal calls:" << addData_call_count
                     << "\nAverage time:" << avg_time * 1000 << "ms"
                     << "\nMinimum total time:" << min_addData_time * 1000 << "ms"
                     << "\nMaximum total time:" << max_addData_time * 1000 << "ms";
        }
    }

    // TTL functions
    int connectTriggerPort_TTL();
    void send_trigger_TTL();
    void set_enable_TTL(bool status);

    void setROIMeans(const Eigen::VectorXd& means);

    Eigen::VectorXd getROIMeanData(int roiIndex);

    void setROINames(const std::vector<std::string>& names) { ROI_names = names; }

signals:
    void channelNamesUpdated(const std::vector<std::string>& channelNames);
    void channelDataUpdated(int channel, const Eigen::VectorXd &data, 
                           const Eigen::VectorXi &triggers_A, 
                           const Eigen::VectorXi &triggers_B, 
                           const Eigen::VectorXi &triggers_out, 
                           const Eigen::VectorXd &time_stamps, 
                                    const size_t &data_index,
                                     std::string source_name);

public slots:
    void savePreprocessingOutput(const Eigen::MatrixXd &output);

private slots:
    void updateSignalViewerData();

private:
    HandlerState handler_state = WAITING_FOR_START;

    QTimer *timer = nullptr;

    // Synchronization primitives
    std::mutex dataMutex;
    std::condition_variable data_condition;
    bool processingWorkerRunning = true;

    int triggerSource;
    Eigen::VectorXi source_channels_;
    bool channel_names_set = false;
    std::vector<std::string> channel_names_;
    Eigen::MatrixXd sample_buffer_;
    Eigen::VectorXd time_stamp_buffer_;
    Eigen::VectorXi trigger_buffer_A;
    Eigen::VectorXi trigger_buffer_B;
    Eigen::VectorXi trigger_buffer_out;
    size_t current_data_index_ = 0;
    int current_sequence_number_ = 0;
    int last_retrieved_sequence_number_ = -1;
    int buffer_length_in_seconds_ = 30;
    int buffer_capacity_;
    int channel_count_;
    int sampling_rate_;
    int simulation_delivery_rate_;

    Eigen::VectorXd processing_sample_vector;

    Eigen::MatrixXd preprocessing_output;
    Eigen::VectorXi preprocessing_triggers_A;
    Eigen::VectorXi preprocessing_triggers_B;
    Eigen::VectorXi preprocessing_triggers_out;
    Eigen::VectorXd preprocessing_time_stamps;
    int preprocessing_number_of_samples;
    int processing_sequence_number;
    int processing_downsampling_factor = 1;
    int prep_buffer_capacity_;

    // Signal viewer data and parameters
    int number_of_EEG_channels = 5;                 //THIS MIGHT NEED SOME INITIALIZATION
    Eigen::MatrixXd PreProcessing_output_save;
    int save_index_tracker = 0;

    // Gradient artifact correction
    bool Apply_GACorr = false;
    GACorrection GACorr_;
    int TA_length = 5000;
    int TR_length = 10000;
    int GA_average_length = 25;
    int TA_tracker = -1;
    bool TA_in_progress = false;
    bool TR_in_progress = false;
    bool GA_is_continuous = false;

    int GA_tracker = -1;
    int GA_shift_front = 100;
    int GA_shift_back = 100;
    bool GA_in_progress = false;

    // Baseline correction
    bool Apply_baseline = false;
    Eigen::VectorXd baseline_average;
    Eigen::VectorXd GA_sum;
    Eigen::VectorXd GA_average;

    // Geometric sum correction
    bool Apply_geometric_sum = false;
    Eigen::VectorXd baseline_correction;
    bool baseline_reset = false;
    double baseline_update_rate = 0.0006;

    // Filter
    bool Apply_filter = false;
    std::vector<double> filterCoeffs_;
    MultiChannelRealTimeFilter RTfilter_;

    // Sending triggers
    TMSConnectionType TMS_connectionType = COM;

    // MAGPRO
    magPro magPro_3G;
    bool triggerPortState = false;
    bool triggerEnableState = false;

    // LABJACK
    int LJM_dtANY = 0;
    int LJM_ctANY = 0;
    int labjack_handle = -1;

    // Time limit in milliseconds
    int time_limit = 500;
    const int min_time_limit = 100;
    const int max_time_limit = 100000;
    std::chrono::time_point<std::chrono::system_clock> latest_trigger_time;

    std::mutex triggerMutex;
    std::unordered_set<int> triggerSet;

    // int SAVE_INDEX_TRACKER = 0;
    // bool data_saved = false;
    // Eigen::MatrixXd save_matrix = Eigen::MatrixXd::Zero(12, 1000000);
    
    int last_save_index = -1;
    bool data_saved = false;
    std::vector<int> seqNum_list;

    Eigen::MatrixXd sample_buffer_save;
    int sample_buffer_save_index = 0;

    Eigen::VectorXd ROI_fMRI_means;
    Eigen::MatrixXd ROI_means_save;  // Add this new matrix to store ROI means over time

    std::vector<std::string> ROI_names;

    // Timing variables for addData
    std::chrono::duration<double> total_addData_time{0};
    std::chrono::duration<double> min_addData_time{std::numeric_limits<double>::max()};
    std::chrono::duration<double> max_addData_time{std::numeric_limits<double>::min()};
    int addData_call_count = 0;
};

#endif // DATAHANDLER_H
