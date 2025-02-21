#ifndef PHASEESTWINDOW_H
#define PHASEESTWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QDockWidget>

#include "processingglwidget.h"
#include "polarHistogramOpenGLWidget.h"
#include "workers/phaseEstimationWorker.h"
#include "../dataHandler/dataHandler.h"

namespace Ui {
class phaseEstwindow;
}

class phaseEstwindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit phaseEstwindow(dataHandler &handler, 
               volatile std::sig_atomic_t &processingWorkerRunning, 
                          Eigen::MatrixXd &processed_data, 
                                  QWidget *parent = nullptr);
    ~phaseEstwindow();

    ProcessingGlWidget* getProcessingGlWidget() { return processingglWidget; }
    PolarHistogramOpenGLWidget* getHistogramWidget() { return HistogramWidget; }

signals:
    void setCustomScaleStatus(bool status);
    void setCustomScaleMin(double min);
    void setCustomScaleMax(double max);

    bool getCustomScaleStatus();
    double getCustomScaleMin();
    double getCustomScaleMax();
    void setFilterState(bool isChecked);
    void setEEGViewState(bool isChecked);
    void setPhaseTargetingState(bool isChecked);
    void setphaseEstimateState(bool isChecked);
    void setSpatilaTargetChannel(int index);

    void setSNRcheck(bool isChecked);
    void sendSNRmax(double value);
    void setSNRthreshold(double value);

    void setPhaseEstParams(phaseEstimateParameters phaseEstParams);

    void setShowTriggers_A(bool isChecked);
    void setShowTriggers_B(bool isChecked);
    void setShowTriggers_out(bool isChecked);
    void setPhaseError(bool isChecked);
    void switchPause();
    void setDrawXaxis(bool isChecked);
    void updateTLineSpacing(int value);
    void updateWidgetChannelNames(std::vector<std::string> processing_channel_names);

    void outerElectrodesStateChanged(std::vector<bool> outerElectrodeCheckStates);
    void setPhaseErrorType(int index);

    void requestEstStates();

public slots:
    // void updateData();
    void updateSpatialChannelNames(std::vector<std::string> names);
    void setNumSamples(int numSamples) { if (processingglWidget) processingglWidget->updateWindowLength_seconds(numSamples / sampling_rate_); }
    void newEstStates(phaseEstimateStates states);
    void newSNRmax(double value);
    void newSNRmax_list(const std::vector<double>& list);
    void showOpenGLDock() { dockWidget->show(); }

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
    void on_checkBox_PhaseTargeting_stateChanged(int arg1);
    void on_checkBox_Stimulation_stateChanged(int arg1);
    void on_checkBox_phaseEstimate_stateChanged(int arg1);
    void on_setParamsButton_clicked();
    void on_comboBox_spatialTarget_currentIndexChanged(int index);
    void setupComboBox_OuterElectrode(int spatial_target_index);
    void handleCheckboxChange(QStandardItem* item);
    void on_pushButton_pause_view_clicked();
    void on_checkBox_triggers_A_stateChanged(int arg1);
    void on_checkBox_triggers_B_stateChanged(int arg1);
    void on_checkBox_phaseError_stateChanged(int arg1);
    void on_checkBox_Xaxis_stateChanged(int arg1);
    void on_lineEdit_XaxisSpacing_editingFinished();
    void on_checkBox_triggers_out_stateChanged(int arg1);
    void on_comboBox_phaseError_currentIndexChanged(int index);

    void on_checkBox_SNRcheck_stateChanged(int arg1);

    void on_SNRtreshold_editingFinished();

    void on_SNRmax_editingFinished();

private:
    Ui::phaseEstwindow *ui;
    ProcessingGlWidget *processingglWidget;
    QDockWidget *dockWidget;
    PolarHistogramOpenGLWidget *HistogramWidget;

    dataHandler &handler;

    phaseEstimateParameters phaseEstParams;
    std::vector<std::string> spatial_channel_names;
    QStringList outer_channel_names;
    Eigen::MatrixXd &processed_data;
    std::mutex dataMutex;

    int numOuterElectrodes = 4;
    std::vector<bool> outerElectrodeCheckStates_;

    int numSamples_ = 10000;
    int sampling_rate_ = 5000;

    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // PHASEESTWINDOW_H
