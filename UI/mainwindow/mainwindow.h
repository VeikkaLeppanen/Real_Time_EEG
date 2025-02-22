#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QToolBar>
#include <QStatusBar>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QVector>
#include <QListWidget>
#include <QListWidgetItem>
#include <QSlider>
#include <QColorDialog>
#include <QLabel>
#include <QMessageBox>
#include <QLineEdit>
#include <QIntValidator>
#include <iostream>
#include "../../workers/preProcessingWorker.h"
#include "../../workers/phaseEstimationWorker.h"
#include "../../workers/EEGSpinWorker.h"
#include "mainglwidget.h"
#include "customTitleBar.h"
#include "../dataHandler/dataHandler.h"
#include "eegwindow/eegwindow.h"
#include "phaseEstimationwindow/phaseEstwindow.h"
#include "TMSwindow/TMSwindow.h"
#include "MRIwindow/mriwindow.h"
#include "MRIwindow/openglMRIwidget.h"


QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

protected:
    bool eventFilter(QObject *obj, QEvent *event) override;
    void closeEvent(QCloseEvent *event) override;

public:
    MainWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent = nullptr);
    ~MainWindow();
    
    void toggleEegBridgeCondition(bool isOn) {
        eegBridgeLabel->setText(
            QString("EEG Bridge: <span style='color: %1;'>%2</span>")
                .arg(isOn ? "green" : "red")
                .arg(isOn ? "ON" : "OFF")
        );
    }

    void togglePreprocessingCondition(bool isOn) {
        preprocessingLabel->setText(
            QString("EEG preprocessing: <span style='color: %1;'>%2</span>")
                .arg(isOn ? "green" : "red")
                .arg(isOn ? "ON" : "OFF")
        );
    }

    void togglePhaseEstCondition(bool isOn) {
        phaseEstLabel->setText(
            QString("EEG phase estimation: <span style='color: %1;'>%2</span>")
                .arg(isOn ? "green" : "red")
                .arg(isOn ? "ON" : "OFF")
        );
    }

    void togglefMRICondition(bool isOn) {
        fMRILabel->setText(
            QString("fMRI: <span style='color: %1;'>%2</span>")
                .arg(isOn ? "green" : "red")
                .arg(isOn ? "ON" : "OFF")
        );
    }

signals:

public slots:
    void updateData();
    void eegBridgeSpin(int port, int timeout);
    void setGACorrection(int GALength, int GAAverage);
    void startGACorrection();
    void stopGACorrection();
    void updateEEGChannels(const std::vector<std::string>& channelNames);
    void updateMRIROIs();

    void startPreprocessing(preprocessingParameters &parameters);
    void startPhaseEstimationprocessing(phaseEstimateParameters phaseEstParams);

    void toggleTMSCondition(bool isOn) {
        TMSLabel->setText(
            QString("TMS: <span style='color: %1;'>%2</span>")
                .arg(isOn ? "green" : "red")
                .arg(isOn ? "ON" : "OFF")
        );
    }

private slots:
    void handleError(const QString& error);

    void EEG_clicked();
    void resetEegWindowPointer();

    void processing_clicked();
    void connect_processing_worker();
    void connect_EEG_Prepworker();
    void resetProcessingWindowPointer();

    void triggering_clicked();
    void resetTMSwinPointer();

    void MRI_clicked();
    void resetMRIwinPointer();

    void addSignalViewer();
    
    void updateSignalViewers();

private:
    Ui::MainWindow *ui;
    QToolBar *mainToolbar;
    QLabel *eegBridgeLabel;
    QLabel *preprocessingLabel;
    QLabel *phaseEstLabel;
    QLabel *fMRILabel;
    QLabel *TMSLabel;

    QVector<QDockWidget*> signalViewers;
    QVector<MainGlWidget*> glWidgets;
    QStringList signalSources = {"(empty)"};

    // eeg_bridge parameters
    EegBridge bridge;

    // Graph parameters
    int samples_to_display = 10000;
    Eigen::MatrixXd processed_data;
    std::vector<std::string> processing_channel_names;

    // Handler parameters
    dataHandler &handler;

    bool preprocessingWorkerRunning = false;
    bool phaseEstWorkerRunning = false;

    // windows
    eegWindow *eegwindow = nullptr;
    phaseEstwindow *phaseEstwin = nullptr;
    TMSwindow *TMSwin = nullptr;
    mriWindow *MRIwin = nullptr;
    preProcessingWorker *preProcessingworker = nullptr;
    phaseEstimationWorker *phaseEstworker = nullptr;

    phaseEstimateParameters phaseEstParams;
    
    volatile std::sig_atomic_t &signal_received;
    volatile std::sig_atomic_t processingWorkerRunning = 0;

    // Testing parameters
    double total_time = 0.0;
    int count = 0;
    Eigen::MatrixXd EEG_output;
};
#endif // MAINWINDOW_H
