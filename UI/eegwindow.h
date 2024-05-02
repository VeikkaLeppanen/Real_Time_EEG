#ifndef EEGWINDOW_H
#define EEGWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QMessageBox>
#include <QObject>
#include <array>
#include "worker.h"
#include "glwidget.h"
#include "eegwindow.h"
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
    explicit eegWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent = nullptr);
    ~eegWindow();

public slots:
    void updateData();

private slots:

    void on_connectButton_clicked();
    void handleError(const QString& error);

    void on_disconnectButton_clicked();

private:
    Ui::EegWindow *ui;

    // eeg_bridge parameters
    EegBridge bridge;
    int port = 50000;

    // Handler parameters
    dataHandler &handler;
    int GALength = 5000;
    int GAAverage = 25;
    int samples_to_display = 30000;

    volatile std::sig_atomic_t &signal_received;
};

#endif // EEGWINDOW_H
