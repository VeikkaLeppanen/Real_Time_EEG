#ifndef PROCESSINGWINDOW_H
#define PROCESSINGWINDOW_H

#include <QMainWindow>
#include <QMessageBox>

#include "processingglwidget.h"
#include "processingworker.h"
#include "../dataHandler/dataHandler.h"

namespace Ui {
class ProcessingWindow;
}

class ProcessingWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit ProcessingWindow(dataHandler &handler, volatile std::sig_atomic_t &processingWorkerRunning, Eigen::MatrixXd &processed_data, QWidget *parent = nullptr);
    ~ProcessingWindow();

signals:
    void startProcessing(processingParameters& params);

public slots:
    void updateData();
    void on_startButton_clicked();
    void on_stopButton_clicked();

private slots:
    void on_downsampling_editingFinished();

    void on_delay_editingFinished();

    void on_edge_editingFinished();

    void on_modelOrder_editingFinished();

    void on_hilbertLength_editingFinished();

    void on_stimulationTarget_editingFinished();

    void on_phaseShift_editingFinished();

private:
    Ui::ProcessingWindow *ui;

    dataHandler &handler;

    processingParameters params;
    Eigen::MatrixXd &processed_data;

    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // PROCESSINGWINDOW_H
