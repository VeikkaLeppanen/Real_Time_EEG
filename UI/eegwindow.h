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
#include <QString>
#include <QDir>
#include <QDebug>

#include <array>

#include "worker.h"
#include "glwidget.h"
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

signals:
    void connectEegBridge(int port);
    void applyGACorrection(int GALength, int GAAverage);
    void stopGACorrection();
    void updateChannelNames(QStringList channelNames);

private slots:

    void on_connectButton_clicked();
    void handleError(const QString& error);

    void on_disconnectButton_clicked();

    void on_lineEditPort_editingFinished();

    void on_lineEditGALength_editingFinished();
    void on_lineEditGAaverage_editingFinished();

    void on_GACorrectionStart_clicked();
    void on_GACorrectionStop_clicked();

    void on_sourceChannelLoad_clicked();
    
private:
    Ui::EegWindow *ui;

    // eeg_bridge parameters
    int port = 50000;

    // Handler parameters
    dataHandler &handler;
    int GALength = 5000;
    int GAAverage = 25;
    int samples_to_display = 30000;

    volatile std::sig_atomic_t &signal_received;
};

#endif // EEGWINDOW_H
