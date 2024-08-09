#include "./mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::MainWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::MainWindow),
      handler(handler),
      signal_received(signal_received)
{
    ui->setupUi(this);
    resize(1280, 720);

    MainGlWidget* mainglWidget = ui->mainGlWidget;
    if (false && mainglWidget) {                                        // REMOVE FALSE IN ORDER TO ENABLE THE GLWIDGET

        connect(mainglWidget, &MainGlWidget::fetchData, this, &MainWindow::updateData);

    } else {
        // Error handling if glWidget is not found
        qWarning("Glwidget not found in UI!");
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateData()
{
    MainGlWidget* mainglWidget = ui->mainGlWidget;
    if (mainglWidget && processingWorkerRunning && (processed_data.size() > 0)) {

        mainglWidget->updateMatrix(processed_data);
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
    handler.GACorr_off();
    handler.reset_GACorr_tracker();
}

void MainWindow::on_EEG_clicked()
{
    if (!eegwindow) {
        eegwindow = new eegWindow(handler, signal_received, processingWorkerRunning, this);
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
    connect(eegwindow, &eegWindow::startPreprocessing, this, &MainWindow::startPreprocessing);
    
    MainGlWidget* mainglWidget = ui->mainGlWidget;
    if (mainglWidget) {
        // connect(eegwindow, &ProcessingWorker::updateChannelNamesQt, mainglWidget, &MainGlWidget::updateChannelNamesQt);
    } else {
        qWarning("Glwidget not found in UI!");
    }
}

void MainWindow::resetEegWindowPointer() {
    eegwindow = nullptr;  // Reset the pointer after the window is destroyed
}

void MainWindow::on_processing_clicked()
{
    if (!processingWindow) {
        processingWindow = new ProcessingWindow(handler, processingWorkerRunning, processed_data, this);
        processingWindow->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
        connect(processingWindow, &ProcessingWindow::destroyed, this, &MainWindow::resetProcessingWindowPointer);
    }
    processingWindow->show();
    processingWindow->raise();
    processingWindow->activateWindow();

    if (worker) {
        QObject::connect(worker, &ProcessingWorker::updatePhaseEstDisplayedData, processingWindow->getProcessingGlWidget(), &ProcessingGlWidget::updateMatrix);
        QObject::connect(worker, &ProcessingWorker::updatePhaseEstwindowNames, processingWindow->getProcessingGlWidget(), &ProcessingGlWidget::updateChannelNamesSTD);
        QObject::connect(worker, &ProcessingWorker::updateSpatialChannelNames, processingWindow, &ProcessingWindow::updateSpatialChannelNames);
        QObject::connect(processingWindow, &ProcessingWindow::setFilterState, worker, &ProcessingWorker::setFilterState);
        QObject::connect(processingWindow, &ProcessingWindow::setEEGViewState, worker, &ProcessingWorker::setEEGViewState);
        QObject::connect(processingWindow, &ProcessingWindow::setphaseEstimateState, worker, &ProcessingWorker::setPhaseEstimationState);
        QObject::connect(processingWindow, &ProcessingWindow::setStimulationState, worker, &ProcessingWorker::setPhaseTargetingState);
        QObject::connect(processingWindow, &ProcessingWindow::setPhaseEstParams, worker, &ProcessingWorker::setPhaseEstimateParameters);
        QObject::connect(processingWindow, &ProcessingWindow::setSpatilaTargetChannel, worker, &ProcessingWorker::setSpatilaTargetChannel);
    }

    // REPLACE WITH A PHASE ESTIMATE BEGIN CONNECT
    // connect(processingWindow, &ProcessingWindow::startProcessing, this, &MainWindow::startPreprocessing);
}

void MainWindow::resetProcessingWindowPointer() {
    processingWindow = nullptr;  // Reset the pointer after the window is destroyed
}

void MainWindow::startPreprocessing(preprocessingParameters& prepParams, phaseEstimateParameters &phaseEstParams)
{
    if (!processingWorkerRunning && handler.isReady()) {
        processingWorkerRunning = 1;

        QThread* thread = new QThread;
        worker = new ProcessingWorker(handler, processed_data, processingWorkerRunning, prepParams, phaseEstParams);
        worker->moveToThread(thread);

        QObject::connect(thread, &QThread::started, worker, &ProcessingWorker::process_start);       // Switch between differentprocessing functions here
        QObject::connect(worker, &ProcessingWorker::finished, thread, &QThread::quit);
        QObject::connect(worker, &ProcessingWorker::error, this, &MainWindow::handleError);
        QObject::connect(worker, &ProcessingWorker::finished, worker, &ProcessingWorker::deleteLater);
        QObject::connect(thread, &QThread::finished, thread, &QThread::deleteLater);

        QObject::connect(thread, &QThread::finished, this, [=]() {
            processingWorkerRunning = 0;
            qDebug("Thread and Worker cleaned up properly");
        });

        thread->start();

        // QObject::connect(worker, &ProcessingWorker::updateProcessingChannelNames, processingWindow, &ProcessingWindow::updateWidgetChannelNames);
        QObject::connect(worker, &ProcessingWorker::updateEEGDisplayedData, eegwindow->getGlWidget(), &Glwidget::updateMatrix);
        QObject::connect(worker, &ProcessingWorker::updateEEGwindowNames, eegwindow->getGlWidget(), &Glwidget::updateChannelNamesSTD);
        QObject::connect(eegwindow, &eegWindow::viewStateChanged, worker, &ProcessingWorker::updateViewState);
        emit eegwindow->viewStateChanged(0);

        if (processingWindow) {
            QObject::connect(worker, &ProcessingWorker::updatePhaseEstDisplayedData, processingWindow->getProcessingGlWidget(), &ProcessingGlWidget::updateMatrix);
            QObject::connect(worker, &ProcessingWorker::updatePhaseEstwindowNames, processingWindow->getProcessingGlWidget(), &ProcessingGlWidget::updateChannelNamesSTD);
            QObject::connect(worker, &ProcessingWorker::updateSpatialChannelNames, processingWindow, &ProcessingWindow::updateSpatialChannelNames);
            QObject::connect(processingWindow, &ProcessingWindow::setFilterState, worker, &ProcessingWorker::setFilterState);
            QObject::connect(processingWindow, &ProcessingWindow::setEEGViewState, worker, &ProcessingWorker::setEEGViewState);
            QObject::connect(processingWindow, &ProcessingWindow::setphaseEstimateState, worker, &ProcessingWorker::setPhaseEstimationState);
            QObject::connect(processingWindow, &ProcessingWindow::setStimulationState, worker, &ProcessingWorker::setPhaseTargetingState);
            QObject::connect(processingWindow, &ProcessingWindow::setPhaseEstParams, worker, &ProcessingWorker::setPhaseEstimateParameters);
            QObject::connect(processingWindow, &ProcessingWindow::setSpatilaTargetChannel, worker, &ProcessingWorker::setSpatilaTargetChannel);
        }

    } else {
        std::cout << "Processing start failed" << '\n';
    }
}

void MainWindow::on_triggering_clicked()
{
    std::cout << "Triggering start" << '\n';
    if (!triggeringWindow) {
        triggeringWindow = new TriggeringWindow(handler, signal_received, this);
        triggeringWindow->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
        connect(triggeringWindow, &TriggeringWindow::destroyed, this, &MainWindow::resetEegWindowPointer);
    }
    triggeringWindow->show();
    triggeringWindow->raise();
    triggeringWindow->activateWindow();
}

void MainWindow::resetTriggeringWindowPointer() {
    triggeringWindow = nullptr;  // Reset the pointer after the window is destroyed
}



