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
#include "../eeg_bridge/eeg_bridge.h"
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

private slots: // Correct usage of slots here
    void on_pushButton_clicked();
    void handleError(const QString& error);
    void on_pushButton_2_clicked();

    void on_lineEditPort_editingFinished();

    void on_lineEditGALength_editingFinished();

    void on_lineEditGAaverage_editingFinished();

    void on_EEG_clicked();
    void resetEegWindowPointer();

private:
    Ui::MainWindow *ui;

    // Handler parameters
    dataHandler &handler;
    int GALength = 5000;
    int GAAverage = 25;

    // eeg_bridge parameters
    EegBridge bridge;
    int port = 50000;

    // eeg window
    eegWindow *eegwindow = nullptr;

    volatile std::sig_atomic_t &signal_received;
};
#endif // MAINWINDOW_H
