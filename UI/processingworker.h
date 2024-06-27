#ifndef PROCESSINGWORKER_H
#define PROCESSINGWORKER_H

#include <QObject>
#include <csignal>
#include <iostream>
#include <cmath>    // For M_PI
#include "../dataHandler/dataHandler.h"
#include "../dataProcessor/processingFunctions.h"


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
    double stimulation_target = M_PI * 0.5;
    int phase_shift = -40;                        // for 5000Hz
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
    explicit ProcessingWorker(dataHandler &handler, Eigen::MatrixXd& processed_data, volatile std::sig_atomic_t &processingWorkerRunning, const processingParameters& params,  QObject* parent = nullptr);
    ~ProcessingWorker();

signals:
    void finished();
    void error(QString err);

public slots:
    void process();
    void process_testing();
    void process_ar_testing();

private:
    dataHandler &handler;
    Eigen::MatrixXd &processed_data;
    volatile std::sig_atomic_t &processingWorkerRunning;
    const processingParameters &params;

    // Filtering parameters
    std::vector<double> filterCoeffs_;
    std::vector<double> b; 
    std::vector<double> a;
};

#endif // PROCESSINGWORKER_H
