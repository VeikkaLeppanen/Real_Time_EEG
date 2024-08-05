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
    void updateProcessingChannelNames(std::vector<std::string> processing_channel_names);
    void updateDisplayedData(const Eigen::MatrixXd &newMatrix);

public slots:
    void process_start() {
        process_future = QtConcurrent::run([this]() { process(); });
    };
    void updateViewState(int index) {
        display_state = static_cast<displayState>(index);
    };

private:
    void process();
    void process_testing();
    void process_stimulation_testing();


    dataHandler &handler;
    Eigen::MatrixXd &processed_data;
    std::mutex dataMutex;
    volatile std::sig_atomic_t &processingWorkerRunning;
    preprocessingParameters &prepParams;
    phaseEstimateParameters &phaseEstParams;
    QFuture<void> process_future;

    bool performPreprocessing = true;
    bool performPhaseEstimation = false;

    displayState display_state = displayState::CWL;

    const bool debug = false;
    void print_debug(std::string msg) {
        if (debug) std::cout << msg << std::endl;
    };
};

#endif // PROCESSINGWORKER_H
