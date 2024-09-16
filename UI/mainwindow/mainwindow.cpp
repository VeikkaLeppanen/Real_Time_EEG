#include "./mainwindow.h"
#include "../ui_mainwindow.h"

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
        thread->setPriority(QThread::TimeCriticalPriority);
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

    if (preProcessingworker) {
        connect_EEG_Prepworker();
        emit eegwindow->requestEstStates();
    }
}

void MainWindow::resetEegWindowPointer() {
    eegwindow = nullptr;  // Reset the pointer after the window is destroyed
}

void MainWindow::startPreprocessing(preprocessingParameters& prepParams)
{
    if (!preprocessingWorkerRunning && handler.isReady()) {
        preprocessingWorkerRunning = true;
        processingWorkerRunning = 1;

        QThread* thread = new QThread;
        preProcessingworker = new preProcessingWorker(handler, processingWorkerRunning, prepParams);
        preProcessingworker->moveToThread(thread);

        QObject::connect(thread, &QThread::started, preProcessingworker, &preProcessingWorker::process_start);       // Switch between differentprocessing functions here
        QObject::connect(preProcessingworker, &preProcessingWorker::finished, thread, &QThread::quit);
        QObject::connect(preProcessingworker, &preProcessingWorker::error, this, &MainWindow::handleError);
        QObject::connect(preProcessingworker, &preProcessingWorker::finished, preProcessingworker, &preProcessingWorker::deleteLater);
        QObject::connect(thread, &QThread::finished, thread, &QThread::deleteLater);

        QObject::connect(thread, &QThread::finished, this, [=]() {
            preProcessingworker = nullptr;
            preprocessingWorkerRunning = false;
            processingWorkerRunning = 0;
            qDebug("Processing Thread and preProcessingWorker cleaned up properly");
        });

        thread->start();
        thread->setPriority(QThread::HighPriority);
        std::cout << "Preprocessing thread start" << '\n';

        connect_EEG_Prepworker();

    } else {
        std::cout << "Preprocessing start failed" << '\n';
    }
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
    
        if (phaseEstworker) connect_processing_worker();
        emit phaseEstwin->requestEstStates();
        
        //Start phase estimation thread
        startPhaseEstimationprocessing(phaseEstParams);
    } else {
        QMessageBox::warning(this, "EEG error", "Please connect the system form EEG windows device tab and start the processing thread from the preprocessing tab.");
    }
}

void MainWindow::connect_processing_worker()
{
    QObject::connect(phaseEstworker, &phaseEstimationWorker::updatePhaseEstDisplayedData, phaseEstwin->getProcessingGlWidget(), &ProcessingGlWidget::updateMatrix);
    QObject::connect(phaseEstworker, &phaseEstimationWorker::updatePhaseEstwindowNames, phaseEstwin->getProcessingGlWidget(), &ProcessingGlWidget::updateChannelNamesSTD);
    QObject::connect(phaseEstworker, &phaseEstimationWorker::updateSpatialChannelNames, phaseEstwin, &phaseEstwindow::updateSpatialChannelNames);
    QObject::connect(phaseEstwin, &phaseEstwindow::setFilterState, phaseEstworker, &phaseEstimationWorker::setFilterState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setEEGViewState, phaseEstworker, &phaseEstimationWorker::setEEGViewState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setphaseEstimateState, phaseEstworker, &phaseEstimationWorker::setPhaseEstimationState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseTargetingState, phaseEstworker, &phaseEstimationWorker::setPhaseTargetingState);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseEstParams, phaseEstworker, &phaseEstimationWorker::setPhaseEstimateParameters);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseError, phaseEstworker, &phaseEstimationWorker::setPhaseDifference);
    QObject::connect(phaseEstwin, &phaseEstwindow::setSpatilaTargetChannel, phaseEstworker, &phaseEstimationWorker::setSpatilaTargetChannel);
    QObject::connect(phaseEstwin, &phaseEstwindow::outerElectrodesStateChanged, phaseEstworker, &phaseEstimationWorker::outerElectrodesStateChanged);
    QObject::connect(phaseEstwin, &phaseEstwindow::setPhaseErrorType, phaseEstworker, &phaseEstimationWorker::setPhaseErrorType);
    QObject::connect(phaseEstworker, &phaseEstimationWorker::sendNumSamples, phaseEstwin, &phaseEstwindow::setNumSamples);
    QObject::connect(phaseEstwin, &phaseEstwindow::requestEstStates, phaseEstworker, &phaseEstimationWorker::sendEstStates);
    QObject::connect(phaseEstworker, &phaseEstimationWorker::newEstStates, phaseEstwin, &phaseEstwindow::newEstStates);
}

void MainWindow::connect_EEG_Prepworker()
{
    QObject::connect(preProcessingworker, &preProcessingWorker::updateEEGDisplayedData, eegwindow->getGlWidget(), &Glwidget::updateMatrix);
    QObject::connect(eegwindow, &eegWindow::setRemoveBCG, preProcessingworker, &preProcessingWorker::setRemoveBCG);
}

void MainWindow::resetProcessingWindowPointer() {
    phaseEstwin = nullptr;  // Reset the pointer after the window is destroyed
}

void MainWindow::startPhaseEstimationprocessing(phaseEstimateParameters phaseEstParams)
{
    if (!phaseEstWorkerRunning && handler.isReady()) {
        phaseEstWorkerRunning = true;
        processingWorkerRunning = 1;

        QThread* thread = new QThread;
        phaseEstworker = new phaseEstimationWorker(handler, processingWorkerRunning, phaseEstParams);
        phaseEstworker->moveToThread(thread);

        QObject::connect(thread, &QThread::started, phaseEstworker, &phaseEstimationWorker::process_start);       // Switch between differentprocessing functions here
        QObject::connect(phaseEstworker, &phaseEstimationWorker::finished, thread, &QThread::quit);
        QObject::connect(phaseEstworker, &phaseEstimationWorker::error, this, &MainWindow::handleError);
        QObject::connect(phaseEstworker, &phaseEstimationWorker::finished, phaseEstworker, &phaseEstimationWorker::deleteLater);
        QObject::connect(thread, &QThread::finished, thread, &QThread::deleteLater);

        QObject::connect(thread, &QThread::finished, this, [=]() {
            phaseEstworker = nullptr;
            phaseEstWorkerRunning = false;
            processingWorkerRunning = 0;
            qDebug("Processing Thread and phaseEstimationWorker cleaned up properly");
        });

        thread->start();
        thread->setPriority(QThread::HighPriority);
        std::cout << "Phase estimation thread start" << '\n';

        if (phaseEstwin) connect_processing_worker();

    } else {
        std::cout << "Phase estimation processing start failed" << '\n';
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



