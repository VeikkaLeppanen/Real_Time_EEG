#ifndef EEGWINDOW_H
#define EEGWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QMessageBox>
#include <QObject>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QIntValidator>
#include <QString>
#include <QDir>
#include <QDebug>
#include <QStandardItemModel>
#include <QStandardItem>
#include <QTimer>

#include <array>

#include "workers/EEGSpinWorker.h"
#include "glwidget.h"
#include "workers/phaseEstimationWorker.h"
#include "../ui_eegwindow.h"
#include "../devices/EEG/eeg_bridge/eeg_bridge.h"
#include "../dataHandler/dataHandler.h"


QT_BEGIN_NAMESPACE
namespace Ui {
class EegWindow;
}
QT_END_NAMESPACE

class eegWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit eegWindow(dataHandler &handler, 
        volatile std::sig_atomic_t &signal_received, 
        volatile std::sig_atomic_t &processingWorkerRunning, 
                           QWidget *parent = nullptr);
    
    ~eegWindow();

    Glwidget* getGlWidget() { return glWidget; }

public slots:
    void updateData();
    void on_startButton_clicked();
    void on_stopButton_clicked();
    void newEstStates(phaseEstimateStates states);

signals:
    void connectEegBridge(int port, int timeout);
    void updateChannelDisplayState(std::vector<bool> channelCheckStates);
    void applyGACorrection(int GALength, int GAAverage);
    void startGACorrection();
    void stopGACorrection();
    void setRemoveBCG(bool isChecked);
    void updateChannelNamesSTD(std::vector<std::string> channelNames);
    void updateChannelNamesQt(QStringList channelNames);
    void scaleDrawStateChanged(bool isChecked);

    void startPreprocessing(preprocessingParameters& prepParams, phaseEstimateParameters &phaseEstParams);
    void viewStateChanged(int index);
    void setShowTriggers_A(bool isChecked);
    void setShowTriggers_B(bool isChecked);
    void switchPause();
    void setDrawXaxis(bool isChecked);
    void updateTLineSpacing(int value);

    void requestEstStates();

    void error(QString err);

private slots:
    void handleError(const QString& error);

    void on_connectButton_clicked();
    void checkHandlerReady();  // Slot to periodically check the handler's readiness
    void on_disconnectButton_clicked();
    void on_sourceChannelLoad_clicked();
    void on_lineEditPort_editingFinished();
    void on_lineEditTimeOut_editingFinished();
    void on_lineEditGraphSamples_editingFinished();
    void on_lineEditGALength_editingFinished();
    void on_lineEditGAaverage_editingFinished();
    void on_checkBox_stateChanged(int arg1);
    void on_filter1_stateChanged(int arg1);
    void on_checkBox_2_stateChanged(int arg1);
    void on_downsampling_editingFinished();
    void on_numberOfSamples_editingFinished();
    void on_delay_editingFinished();
    void on_checkBox_GA_stateChanged(int arg1);
    void on_checkBox_triggers_A_stateChanged(int arg1);
    void on_checkBox_triggers_B_stateChanged(int arg1);
    void on_pushButton_pauseView_clicked();
    void on_checkBox_removeBCG_stateChanged(int arg1);
    void on_checkBox_Xaxis_stateChanged(int arg1);
    void on_lineEdit_XaxisSpacing_editingFinished();
    void updateChannelLength(int value);

private:
    Ui::EegWindow *ui;
    QTimer *checkHandlerTimer;
    Glwidget *glWidget;
    
    bool preprocessingWorkerRunning = false;

    // eeg_bridge parameters
    int port = 50000;
    int bridge_timeout = 10;

    // Graph parameters
    int samples_to_display = 10000;

    // Handler parameters
    dataHandler &handler;
    int GALength = 10000;
    int GAAverage = 25;

    // Preprocessing parameters
    preprocessingParameters prepParams;
    phaseEstimateParameters phaseEstParams;

    void setupChannelNames();
    std::vector<std::string> channelMap_ = {"Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2", "F7", "F8", "T7", "T8", "P7", "P8", "Fz",
                                            "Cz	", "Pz", "Oz", "FC1", "FC2", "CP1", "CP2", "FC5", "FC6", "CP5", "CP6", "TP9", "TP10", "POz", "ECG",
                                            "F1", "F2", "C1", "C2", "P1", "P2", "AF3", "AF4", "FC3", "FC4", "CP3", "CP4", "PO3", "PO4", "F5", "F6",
                                            "C5", "C6", "P5", "P6", "AF7", "AF8", "FT7", "FT8", "TP7", "TP8", "PO7", "PO8", "FT9", "FT10", "Fpz",
                                            "CPz", "CWL1", "CWL2", "CWL3", "CWL4", "CWL5", "CWL6", "CWL7"};
    Eigen::VectorXi source_channels_;
    std::vector<bool> channelCheckStates_;
    QStringList channelNames_;

    volatile std::sig_atomic_t &signal_received;
    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // EEGWINDOW_H
