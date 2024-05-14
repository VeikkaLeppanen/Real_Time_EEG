#include "./mainwindow.h"
#include "./ui_mainwindow.h"

#include "../dataProcessor/processingFunctions.h"

MainWindow::MainWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::MainWindow),
      handler(handler),
      signal_received(signal_received)
{
    ui->setupUi(this);

    MainGlWidget* mainglWidget = ui->mainGlWidget;
    if (mainglWidget) {

        connect(mainglWidget, &MainGlWidget::fetchData, this, &MainWindow::updateData);
        // connect(this, &eegWindow::updateChannelNamesSTD, glWidget, &Glwidget::updateChannelNamesSTD);
        // connect(this, &eegWindow::updateChannelNamesQt, glWidget, &Glwidget::updateChannelNamesQt);
        // connect(this, &eegWindow::updateChannelDisplayState, glWidget, &Glwidget::updateChannelDisplayState);
        // connect(this, &eegWindow::scaleDrawStateChanged, glWidget, &Glwidget::scaleDrawStateChanged);

        // emit updateChannelNamesSTD(handler.getChannelNames());
    } else {
        // Error handling if glWidget is not found
        qWarning("Glwidget not found in UI!");
    }

    // Filtering
    double Fs = 5000;  // Sampling frequency
    double Fc = 120;   // Desired cutoff frequency
    int numTaps = 51;  // Length of the FIR filter
    filterCoeffs_ = designLowPassFilter(numTaps, Fs, Fc);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateData()
{
    MainGlWidget* mainglWidget = ui->mainGlWidget;
    if (mainglWidget && handler.isReady()) {
        
        Eigen::MatrixXd all_channels = handler.returnLatestDataInOrder(samples_to_display);

        // Initializing parameters and matrices for algorithms
        int n_eeg_channels = handler.get_channel_count();
    
        // Filtering
        Eigen::MatrixXd EEG_filtered = Eigen::MatrixXd::Zero(n_eeg_channels, samples_to_display);
    
        EEG_filtered = applyFIRFilterToMatrix(all_channels, filterCoeffs_);

        mainglWidget->updateMatrix(EEG_filtered);
    }
}





// eeg window utilities
void MainWindow::eegBridgeSpin(int port, int timeout)
{
    if (!bridge.isRunning()) {
        signal_received = 0;

        bridge.setPort(port);
        bridge.setTimeout(timeout);

        QThread* thread = new QThread;
        Worker* worker = new Worker(bridge, handler, signal_received);
        worker->moveToThread(thread);

        connect(thread, &QThread::started, worker, &Worker::process);
        connect(worker, &Worker::finished, thread, &QThread::quit);
        connect(worker, &Worker::error, this, &MainWindow::handleError);
        connect(worker, &Worker::finished, worker, &Worker::deleteLater);
        connect(thread, &QThread::finished, thread, &QThread::deleteLater);

        connect(thread, &QThread::finished, this, [=]() {
            qDebug("Thread and Worker cleaned up properly");
        });

        thread->start();
    }
}

void MainWindow::handleError(const QString &error)
{
    // Error handling code here
    QMessageBox::critical(this, "Error", error);
}

void MainWindow::setGACorrection(int GALength, int GAAverage) {
    handler.GACorr_off();

    handler.reset_GACorr(GALength, GAAverage);

    handler.GACorr_on();
}

void MainWindow::startGACorrection() {
    handler.reset_GACorr_tracker();
    handler.GACorr_on();
}

void MainWindow::stopGACorrection() {
    handler.reset_GACorr_tracker();
    handler.GACorr_off();
}

void MainWindow::on_pushButton_2_clicked()
{
    signal_received = 1;
}

void MainWindow::on_EEG_clicked()
{
    if (!eegwindow) {
        eegwindow = new eegWindow(handler, signal_received, this);
        eegwindow->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
        connect(eegwindow, &eegWindow::destroyed, this, &MainWindow::resetEegWindowPointer);
    }
    eegwindow->show();
    eegwindow->raise();
    eegwindow->activateWindow();
    connect(eegwindow, &eegWindow::connectEegBridge, this, &MainWindow::eegBridgeSpin);
    connect(eegwindow, &eegWindow::applyGACorrection, this, &MainWindow::setGACorrection);
    connect(eegwindow, &eegWindow::startGACorrection, this, &MainWindow::startGACorrection);
    connect(eegwindow, &eegWindow::stopGACorrection, this, &MainWindow::stopGACorrection);
}

void MainWindow::resetEegWindowPointer() {
    eegwindow = nullptr;  // Reset the pointer after the window is destroyed
}
