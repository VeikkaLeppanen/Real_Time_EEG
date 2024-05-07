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
    void connectEegBridge(int port, int timeout);
    void updateChannelDisplayState(std::vector<bool> channelCheckStates);
    void applyGACorrection(int GALength, int GAAverage);
    void startGACorrection();
    void stopGACorrection();
    void updateChannelNamesSTD(std::vector<std::string> channelNames);
    void updateChannelNamesQt(QStringList channelNames);

private slots:
    void handleError(const QString& error);

    void on_connectButton_clicked();
    void on_disconnectButton_clicked();
    void on_sourceChannelLoad_clicked();
    void setupComboBox();
    void handleCheckboxChange(QStandardItem* item);

    void on_lineEditPort_editingFinished();
    void on_lineEditTimeOut_editingFinished();

    void on_lineEditGraphSamples_editingFinished();

    void on_lineEditGALength_editingFinished();
    void on_lineEditGAaverage_editingFinished();
    void on_HandlerApplyButton_clicked();
    void on_GACorrectionStart_clicked();
    void on_GACorrectionStop_clicked();

    
private:
    Ui::EegWindow *ui;

    // eeg_bridge parameters
    int port = 50000;
    int bridge_timeout = 10;

    // Graph parameters
    int samples_to_display = 30000;

    // Handler parameters
    dataHandler &handler;
    int GALength = 5000;
    int GAAverage = 25;

    void setupChannelNames();
    std::vector<std::string> channelMap_;
    Eigen::VectorXi source_channels_;
    std::vector<bool> channelCheckStates_;
    QStringList channelNames_;

    volatile std::sig_atomic_t &signal_received;
};

#endif // EEGWINDOW_H
