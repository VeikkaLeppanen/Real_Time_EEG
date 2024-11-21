#include "./mainwindow.h"
#include "../ui_mainwindow.h"

MainWindow::MainWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::MainWindow),
      handler(handler),
      signal_received(signal_received)
{
    ui->setupUi(this);
    setMinimumSize(600, 200); // Ensure a reasonable starting size
    resize(1280, 720);

    mainglWidget = ui->mainGlWidget;

    // ----------------------------
    // Create the menu bar
    // ----------------------------
    QMenuBar *menuBar = this->menuBar();

    // Create the File menu and add actions
    QMenu *eegMenu = menuBar->addMenu("EEG");
    QAction *prepAction = eegMenu->addAction("Preprocessing");
    QAction *phaseEstAction = eegMenu->addAction("Phase estimation");

    // Create the Edit menu and add actions
    QMenu *mriMenu = menuBar->addMenu("MRI");
    QAction *mriSettingsAction = mriMenu->addAction("MRI settings");

    // Create the View menu and add actions
    QMenu *tmsMenu = menuBar->addMenu("TMS");
    QAction *tmsSettingsAction = tmsMenu->addAction("Settings");

    // Connect actions to slots
    connect(prepAction, &QAction::triggered, this, &MainWindow::EEG_clicked);
    connect(phaseEstAction, &QAction::triggered, this, &MainWindow::processing_clicked);
    connect(mriSettingsAction, &QAction::triggered, this, &MainWindow::MRI_clicked);
    connect(tmsSettingsAction, &QAction::triggered, this, &MainWindow::triggering_clicked);

    // ----------------------------
    // Create the status bar
    // ----------------------------
    QStatusBar *statusBar = this->statusBar();

    // Create labels for different conditions
    eegBridgeLabel = new QLabel("EEG Bridge: <span style='color: red;'>OFF</span>");
    preprocessingLabel = new QLabel("EEG preprocessing: <span style='color: red;'>OFF</span>");
    phaseEstLabel = new QLabel("EEG phase estimation: <span style='color: red;'>OFF</span>");
    fMRILabel = new QLabel("fMRI: <span style='color: red;'>OFF</span>");
    TMSLabel = new QLabel("TMS: <span style='color: red;'>OFF</span>");

    // Add labels to the status bar
    statusBar->addWidget(eegBridgeLabel);
    statusBar->addWidget(preprocessingLabel);
    statusBar->addWidget(phaseEstLabel);
    statusBar->addWidget(fMRILabel);
    statusBar->addWidget(TMSLabel);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateData()
{
    MainGlWidget* mainglWidget = ui->mainGlWidget;
    if (mainglWidget && processingWorkerRunning && (processed_data.size() > 0)) {

        // mainglWidget->updateMatrix(processed_data);
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

        toggleEegBridgeCondition(true);
        connect(thread, &QThread::finished, this, [=]() {
            toggleEegBridgeCondition(false);
            qDebug("EEG_bridge Thread and EEGSpinWorker cleaned up properly");
        });
        
        thread->start();
        // thread->setPriority(QThread::HighPriority);
        // thread->setPriority(QThread::TimeCriticalPriority);
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

void MainWindow::EEG_clicked()
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
        // emit eegwindow->requestBCGState();
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

        togglePreprocessingCondition(true);
        QObject::connect(thread, &QThread::finished, this, [=]() {
            preProcessingworker = nullptr;
            preprocessingWorkerRunning = false;
            processingWorkerRunning = 0;
            togglePreprocessingCondition(false);
            qDebug("Processing Thread and preProcessingWorker cleaned up properly");
        });

        thread->start();
        // thread->setPriority(QThread::HighPriority);
        std::cout << "Preprocessing thread start" << '\n';

        connect_EEG_Prepworker();

    } else {
        std::cout << "Preprocessing start failed" << '\n';
    }
}

void MainWindow::processing_clicked()
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
        preprocessingParameters tempParams = preProcessingworker->getPreprocessingParameters();
        phaseEstParams.numberOfSamples = tempParams.numberOfSamples;
        phaseEstParams.downsampling_factor = tempParams.downsampling_factor;
        startPhaseEstimationprocessing(phaseEstParams);
    } else {
        QMessageBox::warning(this, "EEG error", "Please connect the system form EEG preprocessing windows device tab.");
    }
}

void MainWindow::connect_processing_worker()
{
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::updatePhaseEstDisplayedData,    phaseEstwin->getProcessingGlWidget(),       &ProcessingGlWidget::updateMatrix);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::updatePhaseEstwindowNames,      phaseEstwin->getProcessingGlWidget(),       &ProcessingGlWidget::updateChannelNamesSTD);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::updateSpatialChannelNames,      phaseEstwin,                                &phaseEstwindow::updateSpatialChannelNames);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setFilterState,                        phaseEstworker,                             &phaseEstimationWorker::setFilterState);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setEEGViewState,                       phaseEstworker,                             &phaseEstimationWorker::setEEGViewState);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setphaseEstimateState,                 phaseEstworker,                             &phaseEstimationWorker::setPhaseEstimationState);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setPhaseTargetingState,                phaseEstworker,                             &phaseEstimationWorker::setPhaseTargetingState);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setPhaseEstParams,                     phaseEstworker,                             &phaseEstimationWorker::setPhaseEstimateParameters);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setPhaseError,                         phaseEstworker,                             &phaseEstimationWorker::setPhaseDifference);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setSpatilaTargetChannel,               phaseEstworker,                             &phaseEstimationWorker::setSpatilaTargetChannel);
    QObject::connect(phaseEstwin,       &phaseEstwindow::outerElectrodesStateChanged,           phaseEstworker,                             &phaseEstimationWorker::outerElectrodesStateChanged);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setPhaseErrorType,                     phaseEstworker,                             &phaseEstimationWorker::setPhaseErrorType);
    QObject::connect(phaseEstwin,       &phaseEstwindow::requestEstStates,                      phaseEstworker,                             &phaseEstimationWorker::sendEstStates);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setSNRcheck,                           phaseEstworker,                             &phaseEstimationWorker::setSNRcheck);
    QObject::connect(phaseEstwin,       &phaseEstwindow::setSNRthreshold,                       phaseEstworker,                             &phaseEstimationWorker::setSNRthreshold);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::sendNumSamples,                 phaseEstwin,                                &phaseEstwindow::setNumSamples);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::newEstStates,                   phaseEstwin,                                &phaseEstwindow::newEstStates);

    //SNR
    QObject::connect(phaseEstwin,       &phaseEstwindow::sendSNRmax,                            phaseEstworker,                             &phaseEstimationWorker::setSNRmax);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::sendSNRmax,                     phaseEstwin,                                &phaseEstwindow::newSNRmax);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::sendSNRmax_list,                phaseEstwin,                                &phaseEstwindow::newSNRmax_list);
    
    QObject::connect(eegwindow,         &eegWindow::sendPrepStates,                             phaseEstworker,                             &phaseEstimationWorker::receivePrepStates);
    QObject::connect(eegwindow,         &eegWindow::set_processing_pause,                       phaseEstworker,                             &phaseEstimationWorker::set_processing_pause);
}

void MainWindow::connect_EEG_Prepworker()
{
    QObject::connect(preProcessingworker, &preProcessingWorker::updateEEGDisplayedData, eegwindow->getGlWidget(), &Glwidget::updateMatrix);
    QObject::connect(preProcessingworker, &preProcessingWorker::newBCGState, eegwindow, &eegWindow::newBCGState);
    // QObject::connect(eegwindow, &eegWindow::requestBCGState, preProcessingworker, &preProcessingWorker::sendBCGState);
    QObject::connect(eegwindow, &eegWindow::sendPrepStates, preProcessingworker, &preProcessingWorker::setPreprocessingParameters);
    QObject::connect(eegwindow, &eegWindow::setRemoveBCG, preProcessingworker, &preProcessingWorker::setRemoveBCG);
    QObject::connect(eegwindow, &eegWindow::set_processing_pause, preProcessingworker, &preProcessingWorker::set_processing_pause);
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
        
        togglePhaseEstCondition(true);
        QObject::connect(thread, &QThread::finished, this, [=]() {
            phaseEstworker = nullptr;
            phaseEstWorkerRunning = false;
            processingWorkerRunning = 0;
            togglePhaseEstCondition(false);
            qDebug("Processing Thread and phaseEstimationWorker cleaned up properly");
        });

        thread->start();
        // thread->setPriority(QThread::HighPriority);
        std::cout << "Phase estimation thread start" << '\n';

        QObject::connect(preProcessingworker, &preProcessingWorker::preprocessingOutputReady, phaseEstworker, &phaseEstimationWorker::handlePreprocessingOutput);

        if (phaseEstwin) connect_processing_worker();

    } else {
        std::cout << "Phase estimation processing start failed" << '\n';
    }
}

void MainWindow::triggering_clicked()
{
    std::cout << "Triggering start" << '\n';
    if (!TMSwin) {
        TMSwin = new TMSwindow(handler, signal_received, this);
        TMSwin->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
        connect(TMSwin, &TMSwindow::destroyed, this, &MainWindow::resetTMSwinPointer);
        connect(TMSwin, &TMSwindow::TMSStatusChanged, this, &MainWindow::toggleTMSCondition);
    }
    TMSwin->show();
    TMSwin->raise();
    TMSwin->activateWindow();
}

void MainWindow::resetTMSwinPointer() {
    TMSwin = nullptr;  // Reset the pointer after the window is destroyed
}

void MainWindow::MRI_clicked()
{
    std::cout << "MRI start" << '\n';
    if (!MRIwin) {
        MRIwin = new mriWindow(this);
        connect(MRIwin, &TMSwindow::destroyed, this, &MainWindow::resetTMSwinPointer);
    }
    MRIwin->show();
    MRIwin->raise();
    MRIwin->activateWindow();
}

void MainWindow::resetMRIwinPointer() {
    MRIwin = nullptr;  // Reset the pointer after the window is destroyed
}