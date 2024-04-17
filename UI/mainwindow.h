#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QMessageBox>
#include <iostream>
#include "worker.h"
#include "../dataHandler/dataHandler.h"
#include "../eeg_bridge/eeg_bridge.h"

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

private:
    Ui::MainWindow *ui;
    dataHandler& handler;
    EegBridge bridge;
    volatile std::sig_atomic_t &signal_received;
};
#endif // MAINWINDOW_H
