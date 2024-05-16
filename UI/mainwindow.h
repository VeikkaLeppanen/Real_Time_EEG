#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QMessageBox>
#include <QLineEdit>
#include <QIntValidator>
#include <iostream>
#include "processingworker.h"
#include "worker.h"
#include "glwidget.h"
#include "mainglwidget.h"
#include "../dataHandler/dataHandler.h"
#include "eegwindow.h"


QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent = nullptr);
    ~MainWindow();

public slots:
    void updateData();
    void eegBridgeSpin(int port, int timeout);
    void setGACorrection(int GALength, int GAAverage);
    void startGACorrection();
    void stopGACorrection();

private slots:
    void handleError(const QString& error);
    void on_pushButton_2_clicked();

    void on_EEG_clicked();
    void resetEegWindowPointer();

    void on_processingStart_clicked();

    void on_processingStop_clicked();

private:
    Ui::MainWindow *ui;

    // eeg_bridge parameters
    EegBridge bridge;

    // Graph parameters
    int samples_to_display = 10000;
    Eigen::MatrixXd processed_data;

    // Filtering parameters
    std::vector<double> filterCoeffs_;
    std::vector<double> b; 
    std::vector<double> a;

    // Handler parameters
    dataHandler &handler;

    // eeg window
    eegWindow *eegwindow = nullptr;

    volatile std::sig_atomic_t &signal_received;
    volatile std::sig_atomic_t processingWorkerRunning = 0;


    // Testing parameters
    double total_time = 0.0;
    int count = 0;
    Eigen::MatrixXd EEG_output;
};
#endif // MAINWINDOW_H
