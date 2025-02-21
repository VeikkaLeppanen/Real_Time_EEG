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
#include <map>

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
    void preProcessing_start();
    void preProcessing_stop();
    void newBCGState(bool state);

signals:
    void connectEegBridge(int port, int timeout);
    void updateChannelDisplayState(std::vector<bool> channelCheckStates);
    void applyGACorrection(int GALength, int GAAverage);
    void startGACorrection();
    void stopGACorrection();
    void setTRLength(int TRLength);
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

    void requestBCGState();
    void sendPrepStates(preprocessingParameters prepParams);
    void set_processing_pause(bool pause);

    void error(QString err);

private slots:
    void handleError(const QString& error);

    void on_connectButton_clicked();
    void checkHandlerReady();  // Slot to periodically check the handler's readiness
    void on_disconnectButton_clicked();
    void on_sourceChannelLoad_clicked();
    void on_lineEditPort_editingFinished();
    void on_lineEditTimeOut_editingFinished();
    void on_lineEditTALength_editingFinished();
    void on_lineEditTRLength_editingFinished();
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

    void on_checkBox_GeoSum_stateChanged(int arg1);

private:
    Ui::EegWindow *ui;
    QTimer *checkHandlerTimer;
    Glwidget *glWidget;
    
    bool preprocessingWorkerRunning = false;

    // eeg_bridge parameters
    int port = 50000;
    int bridge_timeout = 60;

    // Graph parameters
    int samples_to_display = 10000;

    // Handler parameters
    dataHandler &handler;
    int TALength = 5000;
    int TRLength = 10000;
    int GAAverage = 25;

    // Preprocessing parameters
    preprocessingParameters prepParams;
    phaseEstimateParameters phaseEstParams;

    void setupChannelNames();

    std::map<int, std::string> channelMap_ = { 
        {1, "Fp1"}, {2, "Fp2"}, {3, "F3"}, {4, "F4"}, {5, "C3"}, {6, "C4"}, {7, "P3"}, {8, "P4"}, {9, "O1"}, {10, "O2"},
        {11, "F7"}, {12, "F8"}, {13, "T7"}, {14, "T8"}, {15, "P7"}, {16, "P8"}, {17, "Fz"}, {18, "Cz"}, {19, "Pz"}, {20, "Oz"},
        {21, "FC1"}, {22, "FC2"}, {23, "CP1"}, {24, "CP2"}, {25, "FC5"}, {26, "FC6"}, {27, "CP5"}, {28, "CP6"}, {29, "TP9"}, {30, "TP10"},
        {31, "POz"}, {32, "ECG"},
        {41, "F1"}, {42, "F2"}, {43, "C1"}, {44, "C2"}, {45, "P1"}, {46, "P2"}, {47, "AF3"}, {48, "AF4"}, {49, "FC3"}, {50, "FC4"}, 
        {51, "CP3"}, {52, "CP4"}, {53, "PO3"}, {54, "PO4"}, {55, "F5"}, {56, "F6"}, {57, "C5"}, {58, "C6"}, {59, "P5"}, {60, "P6"}, 
        {61, "AF7"}, {62, "AF8"}, {63, "FT7"}, {64, "FT8"}, {65, "TP7"}, {66, "TP8"}, {67, "PO7"}, {68, "PO8"}, {69, "FT9"}, {70, "FT10"}, 
        {71, "Fpz"}, {72, "CPz"}, {73, "CWL1"}, {74, "CWL2"}, {75, "CWL3"}, {76, "CWL4"}, {77, "CWL5"}, {78, "CWL6"}, {79, "CWL7"}
    };

    Eigen::VectorXi source_channels_;
    std::vector<bool> channelCheckStates_;
    QStringList channelNames_;

    volatile std::sig_atomic_t &signal_received;
    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // EEGWINDOW_H
