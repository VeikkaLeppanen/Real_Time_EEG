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

    ProcessingGlWidget* getProcessingGlWidget() { return processingglWidget; }

signals:
    void setCustomScaleStatus(bool status);
    void setCustomScaleMin(double min);
    void setCustomScaleMax(double max);

    bool getCustomScaleStatus();
    double getCustomScaleMin();
    double getCustomScaleMax();
    void setFilterState(bool isChecked);
    void setEEGViewState(bool isChecked);
    void setStimulationState(bool isChecked);
    void setphaseEstimateState(bool isChecked);
    void setSpatilaTargetChannel(int index);

    void setPhaseEstParams(phaseEstimateParameters phaseEstParams);

    void setShowTriggers_A(bool isChecked);
    void setShowTriggers_B(bool isChecked);
    void switchPause();
    void updateWidgetChannelNames(std::vector<std::string> processing_channel_names);

public slots:
    // void updateData();
    void updateSpatialChannelNames(std::vector<std::string> names);

private slots:
    void on_edge_editingFinished();

    void on_modelOrder_editingFinished();

    void on_hilbertLength_editingFinished();

    void on_stimulationTarget_editingFinished();

    void on_phaseShift_editingFinished();

    void on_scaleMax_editingFinished();

    void on_scaleMin_editingFinished();

    void on_checkBox_stateChanged(int arg1);

    void on_Filter_checkbox_stateChanged(int arg1);

    void on_checkBox_Channels_stateChanged(int arg1);

    void on_checkBox_Stimulation_stateChanged(int arg1);

    void on_checkBox_phaseEstimate_stateChanged(int arg1);

    void on_setParamsButton_clicked();

    void on_comboBox_spatialTarget_currentIndexChanged(int index);

    void on_refreshButton_clicked();

    void on_pushButton_pause_view_clicked();

    void on_checkBox_triggers_A_stateChanged(int arg1);

    void on_checkBox_triggers_B_stateChanged(int arg1);

private:
    Ui::ProcessingWindow *ui;
    ProcessingGlWidget *processingglWidget;

    dataHandler &handler;

    phaseEstimateParameters phaseEstParams;
    std::vector<std::string> spatial_channel_names;
    Eigen::MatrixXd &processed_data;
    std::mutex dataMutex;

    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // PROCESSINGWINDOW_H
