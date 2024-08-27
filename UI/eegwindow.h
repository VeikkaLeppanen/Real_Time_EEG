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

#include "worker.h"
#include "glwidget.h"
#include "processingworker.h"
#include "./ui_eegwindow.h"
#include "../eeg_bridge/eeg_bridge.h"
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


private slots:
    void handleError(const QString& error);

    void on_connectButton_clicked();
    void checkHandlerReady();  // Slot to periodically check the handler's readiness

    void on_disconnectButton_clicked();
    void on_sourceChannelLoad_clicked();
    void setupComboBox();
    void handleCheckboxChange(QStandardItem* item);

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

    void on_comboBox_view_currentIndexChanged(int index);

    void on_checkBox_GA_stateChanged(int arg1);

    void on_checkBox_triggers_A_stateChanged(int arg1);

    void on_checkBox_triggers_B_stateChanged(int arg1);

    void on_pushButton_pauseView_clicked();

    void on_checkBox_removeBCG_stateChanged(int arg1);

    void on_checkBox_Xaxis_stateChanged(int arg1);

    void on_lineEdit_XaxisSpacing_editingFinished();

private:
    Ui::EegWindow *ui;
    QTimer *checkHandlerTimer;
    Glwidget *glWidget;

    // eeg_bridge parameters
    int port = 50000;
    int bridge_timeout = 10;

    // Graph parameters
    int samples_to_display = 30000;

    // Handler parameters
    dataHandler &handler;
    int GALength = 10000;
    int GAAverage = 25;

    // Preprocessing parameters
    preprocessingParameters prepParams;
    phaseEstimateParameters phaseEstParams;

    void setupChannelNames();
    std::vector<std::string> channelMap_;
    Eigen::VectorXi source_channels_;
    std::vector<bool> channelCheckStates_;
    QStringList channelNames_;

    volatile std::sig_atomic_t &signal_received;
    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // EEGWINDOW_H
