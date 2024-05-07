#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QMessageBox>
#include <QLineEdit>
#include <QIntValidator>
#include <iostream>
#include "worker.h"
#include "glwidget.h"
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
    void eegBridgeSpin(int port, int timeout);
    void setGACorrection(int GALength, int GAAverage);
    void startGACorrection();
    void stopGACorrection();

private slots:
    void handleError(const QString& error);
    void on_pushButton_2_clicked();

    void on_EEG_clicked();
    void resetEegWindowPointer();

private:
    Ui::MainWindow *ui;

    // eeg_bridge parameters
    EegBridge bridge;

    // Handler parameters
    dataHandler &handler;

    // eeg window
    eegWindow *eegwindow = nullptr;

    volatile std::sig_atomic_t &signal_received;
};
#endif // MAINWINDOW_H
