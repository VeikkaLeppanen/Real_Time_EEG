#ifndef PREPROCESSINGWORKER_H
#define PREPROCESSINGWORKER_H

#include <QObject>
#include <QMutex>
#include <csignal>
#include <iostream>
#include <cstdlib>
#include <cmath>    // For M_PI
#include "../dataHandler/dataHandler.h"
#include "../EEG/preprocessing/preprocessingFunctions.h"
#include "../EEG/preprocessing/removeBCG.h"
#include "../math/dsp.h"
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

inline std::ostream& operator<<(std::ostream& os, const preprocessingParameters& prepParams) {
    os << "Number of Samples: " << prepParams.numberOfSamples
       << "\nDownsampling Factor: " << prepParams.downsampling_factor
       << "\nDelay: " << prepParams.delay;
    return os;
}

class preProcessingWorker : public QObject {
    Q_OBJECT

public:
    explicit preProcessingWorker(dataHandler &handler, 
               volatile std::sig_atomic_t &processingWorkerRunning, 
                  preprocessingParameters &prepParams_in, 
                                  QObject* parent = nullptr);
    ~preProcessingWorker();

signals:
    void finished();
    void error(QString err);
    void updateEEGDisplayedData(const Eigen::MatrixXd &newMatrix, 
                                const Eigen::VectorXi &triggers_A, 
                                const Eigen::VectorXi &triggers_B, 
                                const Eigen::VectorXd &time_stamps,
                                std::vector<std::string> processing_channel_names);

    void preprocessingOutputReady(const Eigen::MatrixXd &output,
                                  const Eigen::VectorXi &triggers_A,
                                  const Eigen::VectorXi &triggers_B,
                                  const Eigen::VectorXi &triggers_out,
                                  const Eigen::VectorXd &time_stamps,
                                  int number_of_samples,
                                  int seq_num);

public slots:
    void process_start() {
        process_future = QtConcurrent::run([this]() { process(); });
    };

    void setRemoveBCG(bool isChecked) { performRemoveBCG = isChecked; }

    void setPreprocessingParameters(preprocessingParameters newParams);

private:
    void process();

    const bool debug = false;
    void print_debug(std::string msg) {
        if (debug) std::cout << msg << std::endl;
    };

    dataHandler &handler;
    QMutex *mutex;
    volatile std::sig_atomic_t &processingWorkerRunning;
    QFuture<void> process_future;

    bool performRemoveBCG = false;

    int n_EEG_channels_to_use = 5;      
    int n_CWL_channels_to_use = 7;
    int n_channels;
    int samples_to_process;
    int downsampling_factor;
    int downsampled_cols;
    int delay;

    // Memory preallocation for preprocessing matrices
    Eigen::MatrixXd all_channels;
    Eigen::MatrixXd EEG_downsampled;
    Eigen::MatrixXd expCWL;
    Eigen::MatrixXd pinvCWL;
    Eigen::MatrixXd EEG_corrected;

    // Input triggers
    Eigen::VectorXi triggers_A;
    Eigen::VectorXi triggers_B;
    Eigen::VectorXi triggers_out;
    Eigen::VectorXd time_stamps;

    Eigen::MatrixXd EEG_win_data_to_display;
};

#endif // PREPROCESSINGWORKER_H
