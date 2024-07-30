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
    explicit ProcessingWindow(dataHandler &handler, 
               volatile std::sig_atomic_t &processingWorkerRunning, 
                          Eigen::MatrixXd &processed_data, 
                                  QWidget *parent = nullptr);
    ~ProcessingWindow();

signals:
    void startProcessing(processingParameters& params);

    void setCustomScaleStatus(bool status);
    void setCustomScaleMin(double min);
    void setCustomScaleMax(double max);

    bool getCustomScaleStatus();
    double getCustomScaleMin();
    double getCustomScaleMax();

    void updateWidgetChannelNames(std::vector<std::string> processing_channel_names);

public slots:
    void updateData();
    void on_startButton_clicked();
    void on_stopButton_clicked();

private slots:
    void on_downsampling_editingFinished();

    void on_numberOfSamples_editingFinished();

    void on_delay_editingFinished();

    void on_edge_editingFinished();

    void on_modelOrder_editingFinished();

    void on_hilbertLength_editingFinished();

    void on_stimulationTarget_editingFinished();

    void on_phaseShift_editingFinished();

    void on_scaleMax_editingFinished();

    void on_scaleMin_editingFinished();

    void on_checkBox_stateChanged(int arg1);

private:
    Ui::ProcessingWindow *ui;

    dataHandler &handler;

    processingParameters params;
    Eigen::MatrixXd &processed_data;
    std::mutex dataMutex;

    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // PROCESSINGWINDOW_H
