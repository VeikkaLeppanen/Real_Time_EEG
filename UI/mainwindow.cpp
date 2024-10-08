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
        EEGSpinWorker* eeg_worker = new EEGSpinWorker(bridge, handler, signal_received);
        eeg_worker->moveToThread(thread);

        connect(thread, &QThread::started, eeg_worker, &EEGSpinWorker::process);
        connect(eeg_worker, &EEGSpinWorker::finished, thread, &QThread::quit);
        connect(eeg_worker, &EEGSpinWorker::error, this, &MainWindow::handleError);
        connect(eeg_worker, &EEGSpinWorker::finished, eeg_worker, &EEGSpinWorker::deleteLater);
        connect(thread, &QThread::finished, thread, &QThread::deleteLater);

        connect(thread, &QThread::finished, this, [=]() {
            qDebug("EEG_bridge Thread and EEGSpinWorker cleaned up properly");
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

    if (processingworker) {
        connect_EEG_worker();
        emit eegwindow->requestEstStates();
    }
}

void MainWindow::resetEegWindowPointer() {
    eegwindow = nullptr;  // Reset the pointer after the window is destroyed
}

void MainWindow::on_processing_clicked()
{
    if (bridge.isRunning() && handler.isReady() && processingWorkerRunning) {
        if (!phaseEstwin) {
            phaseEstwin = new phaseEstwindow(handler, processingWorkerRunning, processed_data, this);
            phaseEstwin->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
            connect(phaseEstwin, &phaseEstwindow::destroyed, this, &MainWindow::resetProcessingWindowPointer);
        }
        phaseEstwin->show();
        phaseEstwin->raise();
        phaseEstwin->activateWindow();
    
        if (processingworker) connect_processing_worker();
        emit phaseEstwin->requestEstStates();
    } else {
        QMessageBox::warning(this, "EEG error", "Please connect the system form EEG windows device tab and start the processing thread from the preprocessing tab.");
    }
}

void MainWindow::connect_processing_worker()
{
    QObject::connect(processingworker, &ProcessingWorker::updatePhaseEstDisplayedData, phaseEstwin->getProcessingGlWidget(), &ProcessingGlWidget::updateMatrix);
    QObject::connect(processingworker, &ProcessingWorker::updatePhaseEstwindowNames, phaseEstwin->getProcessingGlWidget(), &ProcessingGlWidget::updateChannelNamesSTD);
    QObject::connect(processingworker, &ProcessingWorker::updateSpatialChannelNames, phaseEstwin, &phaseEstwindow::updateSpatialChannelNames);
    QObject::connect(phaseEstwin, &phaseEstwindow::setFilterState, processingworker, &ProcessingWorker::setFilterState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setEEGViewState, processingworker, &ProcessingWorker::setEEGViewState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setphaseEstimateState, processingworker, &ProcessingWorker::setPhaseEstimationState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseTargetingState, processingworker, &ProcessingWorker::setPhaseTargetingState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseEstParams, processingworker, &ProcessingWorker::setPhaseEstimateParameters);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseError, processingworker, &ProcessingWorker::setPhaseDifference);
    QObject::connect(phaseEstwin, &phaseEstwindow::setSpatilaTargetChannel, processingworker, &ProcessingWorker::setSpatilaTargetChannel);
    QObject::connect(phaseEstwin, &phaseEstwindow::outerElectrodesStateChanged, processingworker, &ProcessingWorker::outerElectrodesStateChanged);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseErrorType, processingworker, &ProcessingWorker::setPhaseErrorType);
    QObject::connect(processingworker, &ProcessingWorker::sendNumSamples, phaseEstwin, &phaseEstwindow::setNumSamples);
    QObject::connect(phaseEstwin, &phaseEstwindow::requestEstStates, processingworker, &ProcessingWorker::sendEstStates);
    QObject::connect(processingworker, &ProcessingWorker::newEstStates, phaseEstwin, &phaseEstwindow::newEstStates);
}

void MainWindow::connect_EEG_worker()
{
    QObject::connect(processingworker, &ProcessingWorker::updateEEGDisplayedData, eegwindow->getGlWidget(), &Glwidget::updateMatrix);
    QObject::connect(eegwindow, &eegWindow::setRemoveBCG, processingworker, &ProcessingWorker::setRemoveBCG);
    QObject::connect(eegwindow, &eegWindow::requestEstStates, processingworker, &ProcessingWorker::sendEstStates);
    QObject::connect(processingworker, &ProcessingWorker::newEstStates, eegwindow, &eegWindow::newEstStates);
}

void MainWindow::resetProcessingWindowPointer() {
    phaseEstwin = nullptr;  // Reset the pointer after the window is destroyed
}

void MainWindow::startPreprocessing(preprocessingParameters& prepParams, phaseEstimateParameters &phaseEstParams)
{
    if (!processingWorkerRunning && handler.isReady()) {
        processingWorkerRunning = 1;

        QThread* thread = new QThread;
        processingworker = new ProcessingWorker(handler, processed_data, processingWorkerRunning, prepParams, phaseEstParams);
        processingworker->moveToThread(thread);

        QObject::connect(thread, &QThread::started, processingworker, &ProcessingWorker::process_start);       // Switch between differentprocessing functions here
        QObject::connect(processingworker, &ProcessingWorker::finished, thread, &QThread::quit);
        QObject::connect(processingworker, &ProcessingWorker::error, this, &MainWindow::handleError);
        QObject::connect(processingworker, &ProcessingWorker::finished, processingworker, &ProcessingWorker::deleteLater);
        QObject::connect(thread, &QThread::finished, thread, &QThread::deleteLater);

        QObject::connect(thread, &QThread::finished, this, [=]() {
            processingworker = nullptr;
            processingWorkerRunning = 0;
            qDebug("Processing Thread and ProcessingWorker cleaned up properly");
        });

        thread->start();

        connect_EEG_worker();

        if (phaseEstwin) connect_processing_worker();

    } else {
        std::cout << "Processing start failed" << '\n';
    }
}

void MainWindow::on_triggering_clicked()
{
    std::cout << "Triggering start" << '\n';
    if (!TMSwin) {
        TMSwin = new TMSwindow(handler, signal_received, this);
        TMSwin->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
        connect(TMSwin, &TMSwindow::destroyed, this, &MainWindow::resetTMSwinPointer);
    }
    TMSwin->show();
    TMSwin->raise();
    TMSwin->activateWindow();
}

void MainWindow::resetTMSwinPointer() {
    TMSwin = nullptr;  // Reset the pointer after the window is destroyed
}



