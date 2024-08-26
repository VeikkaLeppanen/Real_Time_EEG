#ifndef PROCESSINGWORKER_H
#define PROCESSINGWORKER_H

#include <QObject>
#include <csignal>
#include <iostream>
#include <cstdlib>
#include <cmath>    // For M_PI
#include "../dataHandler/dataHandler.h"
#include "../dataProcessor/processingFunctions.h"
#include <boost/stacktrace.hpp>
#include <QtConcurrent/QtConcurrent>
#include <QFuture>

struct preprocessingParameters {

    // Number of samples to use for the processing
    int numberOfSamples = 10000;

    //  downsampling
    int downsampling_factor = 10;

    // removeBCG
    int delay = 11;
};

struct phaseEstimateParameters {

    // Filter 9-13Hz
    int filter2_length = 250;

    // phase estimate
    size_t edge = 35;
    size_t modelOrder = 15;
    size_t hilbertWinLength = 64;

    // stimulation
    double stimulation_target = 0;          //M_PI * 0.5;
    int phase_shift = 0;                    // for 5000Hz
};

enum class displayState {
    RAW = 0,
    DOWNSAMPLED = 1,
    CWL = 2
};

inline std::ostream& operator<<(std::ostream& os, const preprocessingParameters& prepParams) {
    os << "Number of Samples: " << prepParams.numberOfSamples
       << "\nDownsampling Factor: " << prepParams.downsampling_factor
       << "\nDelay: " << prepParams.delay;
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const phaseEstimateParameters& phaseEstParams) {
    os << "\nEdge: " << phaseEstParams.edge
       << "\nModel Order: " << phaseEstParams.modelOrder
       << "\nHilbert Window Length: " << phaseEstParams.hilbertWinLength
       << "\nStimulation Target: " << phaseEstParams.stimulation_target
       << "\nPhase Shift: " << phaseEstParams.phase_shift;
    return os;
}

class ProcessingWorker : public QObject {
    Q_OBJECT

public:
    explicit ProcessingWorker(dataHandler &handler, 
                          Eigen::MatrixXd &processed_data, 
               volatile std::sig_atomic_t &processingWorkerRunning, 
                  preprocessingParameters &prepParams, 
                  phaseEstimateParameters &phaseEstParams, 
                                  QObject* parent = nullptr);
    ~ProcessingWorker();

signals:
    void finished();
    void error(QString err);
    void updateEEGwindowNames(std::vector<std::string> processing_channel_names);
    void updatePhaseEstwindowNames(std::vector<std::string> processing_channel_names);
    void updateEEGDisplayedData(const Eigen::MatrixXd &newMatrix, const Eigen::VectorXi triggers_A, const Eigen::VectorXi triggers_B);
    void updatePhaseEstDisplayedData(const Eigen::MatrixXd &newMatrix, const Eigen::VectorXi triggers_A, const Eigen::VectorXi triggers_B, int numPastElements, int numFutureElements);

    void updateSpatialChannelNames(std::vector<std::string> processing_channel_names);

public slots:
    void process_start() {
        process_future = QtConcurrent::run([this]() { process(); });
    };
    void updateViewState(int index) {
        display_state = static_cast<displayState>(index);
    };

    void setPhaseEstimationState(bool isChecked) { performPhaseEstimation = isChecked; }

    void setFilterState(bool isChecked) { performFiltering = isChecked; }
    void setEstimationState(bool isChecked) { performEstimation = isChecked; }
    void setHilbertTransformState(bool isChecked) { performHilbertTransform = isChecked; }
    void setPhaseTargetingState(bool isChecked) { performPhaseTargeting = isChecked; }
    void setEEGViewState(bool isChecked) { phasEst_display_all_EEG_channels = isChecked; }
    void setSpatilaTargetChannel(int index) { spatial_channel_index = index; }

    void outerElectrodesStateChanged(std::vector<bool> outerElectrodeCheckStates) { outerElectrodeCheckStates_ = outerElectrodeCheckStates; };

    void setPhaseEstimateParameters(phaseEstimateParameters newParams) {
        phaseEstParams.edge = newParams.edge;
        phaseEstParams.modelOrder = newParams.modelOrder;
        phaseEstParams.hilbertWinLength = newParams.hilbertWinLength; 
        phaseEstParams.stimulation_target = newParams.stimulation_target;
        phaseEstParams.phase_shift = newParams.phase_shift;

        edge_cut_cols = filter2_length - 2 * newParams.edge;
        estimationLength = newParams.edge + std::ceil(newParams.hilbertWinLength / 2);
        display_length = downsampled_cols + estimationLength - newParams.edge;

        EEG_filter2 = Eigen::VectorXd::Zero(filter2_length);
        EEG_predicted.resize(estimationLength, 0.0);
        EEG_hilbert.resize(estimationLength, std::complex<double>(0.0, 0.0));
        phaseAngles = Eigen::VectorXd::Zero(estimationLength);

        phase_diff_hilbert.resize(estimationLength, std::complex<double>(0.0, 0.0));
        phaseDifference = Eigen::VectorXd::Zero(downsampled_cols - newParams.edge);

        outerElectrodeCheckStates_.resize(numOuterElectrodes, true);

        Data_to_display = Eigen::MatrixXd::Zero(9, display_length);
    }

    Eigen::VectorXd getPhaseDifference_vector() {
        if (phaseDifference_current_index == 0) return phaseDifference;
        
        std::lock_guard<std::mutex> lock(this->dataMutex); // Protect shared data access

        int N = downsampled_cols - phaseEstParams.edge;
        Eigen::VectorXd output(N); // Corrected initialization

        // Calculate the number of samples in each segment
        int right = N - phaseDifference_current_index;

        output.head(right) = phaseDifference.tail(right);
        
        output.tail(phaseDifference_current_index) = phaseDifference.head(phaseDifference_current_index);

        return output;
    }

private:
    void process();


    dataHandler &handler;
    Eigen::MatrixXd &processed_data;
    std::mutex dataMutex;
    volatile std::sig_atomic_t &processingWorkerRunning;
    preprocessingParameters &prepParams;
    phaseEstimateParameters &phaseEstParams;
    QFuture<void> process_future;

    int downsampled_cols;

    bool performPreprocessing = true;
    bool performPhaseEstimation = false;
    bool performSNRcheck = false;
    bool performFiltering = false;
    bool performEstimation = true;
    bool performHilbertTransform = true;
    bool performPhaseTargeting = false;
    bool performPhaseDifference = true;

    int spatial_channel_index = 0;
    int numOuterElectrodes = 4;
    std::vector<bool> outerElectrodeCheckStates_;
    int filter2_length = 250;
    int edge_cut_cols;
    size_t estimationLength;

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

    displayState display_state = displayState::CWL;
    bool phasEst_display_all_EEG_channels = false;
    int display_length;
    Eigen::MatrixXd Data_to_display;

    const bool debug = false;
    void print_debug(std::string msg) {
        if (debug) std::cout << msg << std::endl;
    };
};

#endif // PROCESSINGWORKER_H
