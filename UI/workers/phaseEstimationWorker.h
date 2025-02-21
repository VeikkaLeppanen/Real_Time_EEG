#ifndef PHASEESTIMATIONWORKER_H
#define PHASEESTIMATIONWORKER_H

#include <QObject>
#include <csignal>
#include <iostream>
#include <cstdlib>
#include <cmath>    // For M_PI
#include "../dataHandler/dataHandler.h"
#include "../EEG/preprocessing/preprocessingFunctions.h"
#include "../EEG/preprocessing/removeBCG.h"
#include "../EEG/phaseEstimation/phaseEstimationFunctions.h"
#include "../math/dsp.h"
#include "preProcessingWorker.h"
#include <boost/stacktrace.hpp>
#include <QtConcurrent/QtConcurrent>
#include <QFuture>
#include <chrono>

struct phaseEstimateParameters {

    // Number of samples to use for the processing
    int numberOfSamples = 10000;

    //  downsampling
    int downsampling_factor = 10;

    // Filter 9-13Hz
    int filter2_length = 250;
    double SNR_threshold = 0.3;

    // phase estimate
    size_t edge = 35;
    size_t modelOrder = 15;
    size_t hilbertWinLength = 64;

    // stimulation
    double stimulation_target = 0;          //M_PI * 0.5;    [0, 2*pi]
    int phase_shift = -40;                    // for 5000Hz
};

struct phaseEstimateStates {
    bool performPreprocessing = true;
    bool performPhaseEstimation = false;
    bool performSNRcheck = false;
    bool performRemoveBCG = false;
    bool performFiltering = false;
    bool performEstimation = true;
    bool performHilbertTransform = true;
    bool performPhaseTargeting = false;
    bool performPhaseDifference = false;
    bool phasEst_display_all_EEG_channels = false;
};

inline std::ostream& operator<<(std::ostream& os, const phaseEstimateParameters& phaseEstParams) {
    os << "Number of Samples: " << phaseEstParams.numberOfSamples
       << "\nDownsampling Factor: " << phaseEstParams.downsampling_factor
       << "\nEdge: " << phaseEstParams.edge
       << "\nModel Order: " << phaseEstParams.modelOrder
       << "\nHilbert Window Length: " << phaseEstParams.hilbertWinLength
       << "\nStimulation Target: " << phaseEstParams.stimulation_target
       << "\nPhase Shift: " << phaseEstParams.phase_shift;
    return os;
}

class phaseEstimationWorker : public QObject {
    Q_OBJECT

public:
    explicit phaseEstimationWorker(dataHandler &handler, 
               volatile std::sig_atomic_t &processingWorkerRunning, 
                  phaseEstimateParameters &phaseEstParams_in, 
                                  QObject* parent = nullptr);
    ~phaseEstimationWorker();

signals:
    void finished();
    void error(QString err);

    void updatePhaseEstDisplayedData(const Eigen::MatrixXd &newMatrix, 
                                     const Eigen::VectorXi &triggers_A, 
                                     const Eigen::VectorXi &triggers_B, 
                                     const Eigen::VectorXi &triggers_out, 
                                     const Eigen::VectorXd &time_stamps, 
                                     int numPastElements, 
                                     int numFutureElements);

    void updatePhaseEstwindowNames(std::vector<std::string> processing_channel_names);

    void updateSpatialChannelNames(std::vector<std::string> processing_channel_names);
    void sendNumSamples(int numSamples);
    void newEstStates(phaseEstimateStates states);
    
    void sendSNRmax(double value);
    void sendSNRmax_list(const std::vector<double>& list);

    void polarHistogramAddSample_1(double angle);
    void polarHistogramAddSample_2(double angle);

public slots:
    void process_start() {
        process_future = QtConcurrent::run([this]() { process(); });
    };

    void handlePreprocessingOutput(const Eigen::MatrixXd &output,
                                   const Eigen::VectorXi &triggers_A_in,
                                   const Eigen::VectorXi &triggers_B_in,
                                   const Eigen::VectorXi &triggers_out_in,
                                   const Eigen::VectorXd &time_stamps_in,
                                   int number_of_samples,
                                   int seq_num);

    void setPhaseEstimationState(bool isChecked) { phaseEstStates.performPhaseEstimation = isChecked; }

    void setFilterState(bool isChecked) { phaseEstStates.performFiltering = isChecked; }
    void setEstimationState(bool isChecked) { phaseEstStates.performEstimation = isChecked; }
    void setHilbertTransformState(bool isChecked) { phaseEstStates.performHilbertTransform = isChecked; }
    void setPhaseTargetingState(bool isChecked) { phaseEstStates.performPhaseTargeting = isChecked; }
    void setEEGViewState(bool isChecked) { phaseEstStates.phasEst_display_all_EEG_channels = isChecked; }
    void setPhaseDifference(bool isChecked) { phaseEstStates.performPhaseDifference = isChecked; };
    void setSpatilaTargetChannel(int index) { spatial_channel_index = index; }

    // TODO: SNR check
    void setSNRcheck(bool isChecked) {};
    void setSNRmax(double value) {};
    void setSNRthreshold(double value) {};

    void outerElectrodesStateChanged(std::vector<bool> outerElectrodeCheckStates) { outerElectrodeCheckStates_ = outerElectrodeCheckStates; };

    void setPhaseEstimateParameters(phaseEstimateParameters newParams);

    Eigen::VectorXd getPhaseDifference_vector() {
        if (phaseDifference_current_index == 0) return phaseDifference;
        
        std::lock_guard<std::mutex> lock(this->dataMutex); // Protect shared data access

        int N = downsampled_cols - edge;
        Eigen::VectorXd output(N); // Corrected initialization

        // Calculate the number of samples in each segment
        int right = N - phaseDifference_current_index;

        output.head(right) = phaseDifference.tail(right);
        
        output.tail(phaseDifference_current_index) = phaseDifference.head(phaseDifference_current_index);

        return output;
    }
    void setPhaseErrorType(int index) { phaseErrorType = index; };

    void sendEstStates() { emit newEstStates(phaseEstStates); }
    void receivePrepStates(preprocessingParameters prepParams) {
        currentPhaseEstParams.numberOfSamples = prepParams.numberOfSamples;
        currentPhaseEstParams.downsampling_factor = prepParams.downsampling_factor;
        setPhaseEstimateParameters(currentPhaseEstParams);
    }

    void set_processing_pause(bool pause) { processing_pause = pause; }

private:
    void process();

    const bool debug = false;
    void print_debug(std::string msg) {
        if (debug) std::cout << msg << std::endl;
    };

    dataHandler &handler;
    std::mutex dataMutex;
    volatile std::sig_atomic_t &processingWorkerRunning;
    QFuture<void> process_future;

    int downsampled_cols;

    bool processing_pause = false;
    phaseEstimateParameters currentPhaseEstParams;
    phaseEstimateStates phaseEstStates;

    int spatial_channel_index = 0;
    int numOuterElectrodes = 4;
    std::vector<bool> outerElectrodeCheckStates_;

    int sequence_number = 0;
    int n_EEG_channels_to_use = 5;      
    int n_CWL_channels_to_use = 7;
    int n_channels;
    int samples_to_process;
    int downsampling_factor;
    int delay;
    
    // TODO: SNR check

    size_t edge;
    int filter2_length;
    int edge_cut_cols;
    int estimationLength;
    size_t modelOrder;
    size_t hilbertWinLength; 
    double stimulation_target;
    int phase_shift;

    // Memory preallocation for preprocessing matrices
    Eigen::MatrixXd EEG_corrected;
    Eigen::VectorXd EEG_spatial;

    // Input triggers
    Eigen::VectorXi triggers_A;
    Eigen::VectorXi triggers_B;
    Eigen::VectorXi triggers_out;
    Eigen::VectorXd time_stamps;

    Eigen::VectorXd EEG_filter2;
    std::vector<double> EEG_predicted;
    std::vector<std::complex<double>> EEG_hilbert;
    Eigen::VectorXd phaseAngles;

    // Phase difference/error
    std::vector<std::complex<double>> phase_diff_hilbert;
    Eigen::VectorXd phaseDifference;
    int phaseDifference_current_index = 0;
    double last_phase = 0.0;
    int last_phase_seqnum = -1;

    int phaseErrorType = 0;
    int display_length;
    Eigen::MatrixXd EEG_win_data_to_display;
    Eigen::MatrixXd Data_to_display;

    Eigen::VectorXd Pxx = Eigen::VectorXd::Zero(256);
    Eigen::MatrixXd Pxx_save = Eigen::MatrixXd::Zero(100, 256);
    int Pxx_save_index = 0;
    bool Pxx_save_done = false;

    // Timing variables
    std::chrono::duration<double> total_process_time{0};
    std::chrono::duration<double> min_process_time{std::numeric_limits<double>::max()};
    std::chrono::duration<double> max_process_time{std::numeric_limits<double>::min()};
    int process_call_count = 0;

    // Individual step timing
    std::chrono::duration<double> total_spatial_time{0};
    std::chrono::duration<double> total_filter_time{0};
    std::chrono::duration<double> total_phase_est_time{0};
    int step_call_count = 0;
};

#endif // PHASEESTIMATIONWORKER_H
