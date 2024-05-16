#ifndef PROCESSINGWORKER_H
#define PROCESSINGWORKER_H

#include <QObject>
#include <csignal>
#include "../dataHandler/dataHandler.h"
#include "../dataProcessor/processingFunctions.h"

class dataHandler;

class ProcessingWorker : public QObject {
    Q_OBJECT

public:
    explicit ProcessingWorker(dataHandler &handler, Eigen::MatrixXd& processed_data, volatile std::sig_atomic_t &processingWorkerRunning, QObject* parent = nullptr);
    ~ProcessingWorker();

signals:
    void finished();
    void error(QString err);

public slots:
    void process();

private:
    dataHandler &handler;
    Eigen::MatrixXd &processed_data;
    volatile std::sig_atomic_t &processingWorkerRunning;

    // Filtering parameters
    std::vector<double> filterCoeffs_;
    std::vector<double> b; 
    std::vector<double> a;
    int samples_to_display = 10000;
};

#endif // PROCESSINGWORKER_H
