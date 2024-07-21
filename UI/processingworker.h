#ifndef PROCESSINGWORKER_H
#define PROCESSINGWORKER_H

#include <QObject>
#include <csignal>
#include <iostream>
#include <cmath>    // For M_PI
#include "../dataHandler/dataHandler.h"
#include "../dataProcessor/processingFunctions.h"
#include <boost/stacktrace.hpp>


struct processingParameters {

    // Number of samples to use for the processing
    int numberOfSamples = 10000;

    //  downsampling
    int downsampling_factor = 10;

    // removeBCG
    int delay = 11;

    // phase estimate
    size_t edge = 35;
    size_t modelOrder = 15;
    size_t hilbertWinLength = 64;

    // stimulation
    double stimulation_target = 0; //M_PI * 0.5;
    int phase_shift = 0;                        // for 5000Hz
};

inline std::ostream& operator<<(std::ostream& os, const processingParameters& params) {
    os << "Number of Samples: " << params.numberOfSamples
       << "\nDownsampling Factor: " << params.downsampling_factor
       << "\nDelay: " << params.delay
       << "\nEdge: " << params.edge
       << "\nModel Order: " << params.modelOrder
       << "\nHilbert Window Length: " << params.hilbertWinLength
       << "\nStimulation Target: " << params.stimulation_target
       << "\nPhase Shift: " << params.phase_shift;
    return os;
}

class ProcessingWorker : public QObject {
    Q_OBJECT

public:
    explicit ProcessingWorker(dataHandler &handler, 
                          Eigen::MatrixXd &processed_data, 
               volatile std::sig_atomic_t &processingWorkerRunning, 
               const processingParameters &params, 
                                  QObject* parent = nullptr);
    ~ProcessingWorker();

signals:
    void finished();
    void error(QString err);
    void updateProcessingChannelNames(std::vector<std::string> processing_channel_names);

public slots:
    void process();
    void process_timing();
    void process_testing();
    void process_stimulation_testing();

private:
    dataHandler &handler;
    Eigen::MatrixXd &processed_data;
    volatile std::sig_atomic_t &processingWorkerRunning;
    const processingParameters &params;

};

#endif // PROCESSINGWORKER_H
