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
    // Create the toolbar
    // ----------------------------
    QToolBar *toolbar = addToolBar("Main Toolbar");

    // Create action for adding signal viewer
    QAction *addSignalViewerAction = new QAction("Add Signal Viewer", this);
    toolbar->addAction(addSignalViewerAction);

    // Connect action to slot
    connect(addSignalViewerAction, &QAction::triggered, this, &MainWindow::addSignalViewer);

    setDockNestingEnabled(true);

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

    // Add this line to make signalSources accessible as a property
    setProperty("signalSources", QVariant::fromValue(signalSources));

    // Add this connection
    connect(&handler, &dataHandler::channelNamesUpdated, this, &MainWindow::updateEEGChannels);

    // Add connection for MRI ROI updates
    if (MRIwin && MRIwin->getGlWidget()) {
        connect(MRIwin->getGlWidget(), &OpenGLMRIWidget::ROIListChanged, 
                this, &MainWindow::updateMRIROIs);
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    // Signal all threads to stop
    signal_received = 1;
    
    // Stop EEG Bridge if running
    if (bridge.isRunning()) {
        bridge.running = false;
        bridge.close_socket();
    }

    // Stop preprocessing if running
    if (preProcessingworker) {
        preprocessingWorkerRunning = false;
        processingWorkerRunning = 0;
    }

    // Close all child windows
    if (eegwindow) {
        eegwindow->close();
    }
    if (TMSwin) {
        TMSwin->close();
    }
    if (MRIwin) {
        MRIwin->close();
    }
    if (phaseEstwin) {
        phaseEstwin->close();
    }

    // Clean up signal viewers
    for (auto viewer : signalViewers) {
        viewer->close();
    }
    signalViewers.clear();

    // Give threads time to clean up
    QThread::msleep(100);

    // Accept the close event
    event->accept();
}

void MainWindow::updateData()
{
    // MainGlWidget* mainglWidget = ui->mainGlWidget;
    // if (mainglWidget && processingWorkerRunning && (processed_data.size() > 0)) {

    //     mainglWidget->updateMatrix(processed_data);
    // }
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
            handler.save_seqnum_list();
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
        // eegwindow->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
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
            // phaseEstwin->setAttribute(Qt::WA_DeleteOnClose); // Window is deleted on close
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

    // SNR
    QObject::connect(phaseEstwin,       &phaseEstwindow::sendSNRmax,                            phaseEstworker,                             &phaseEstimationWorker::setSNRmax);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::sendSNRmax,                     phaseEstwin,                                &phaseEstwindow::newSNRmax);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::sendSNRmax_list,                phaseEstwin,                                &phaseEstwindow::newSNRmax_list);
    
    QObject::connect(eegwindow,         &eegWindow::sendPrepStates,                             phaseEstworker,                             &phaseEstimationWorker::receivePrepStates);
    QObject::connect(eegwindow,         &eegWindow::set_processing_pause,                       phaseEstworker,                             &phaseEstimationWorker::set_processing_pause);
    
    // Polar histogram
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::polarHistogramAddSample_1,      phaseEstwin->getHistogramWidget(),          &PolarHistogramOpenGLWidget::addSampleToFirstCircle);
    QObject::connect(phaseEstworker,    &phaseEstimationWorker::polarHistogramAddSample_2,      phaseEstwin->getHistogramWidget(),          &PolarHistogramOpenGLWidget::addSampleToSecondCircle);
}

void MainWindow::connect_EEG_Prepworker()
{
    QObject::connect(preProcessingworker, &preProcessingWorker::updateEEGDisplayedData, eegwindow->getGlWidget(), &Glwidget::updateMatrix);
    QObject::connect(preProcessingworker, &preProcessingWorker::newBCGState, eegwindow, &eegWindow::newBCGState);
    // QObject::connect(eegwindow, &eegWindow::requestBCGState, preProcessingworker, &preProcessingWorker::sendBCGState);
    QObject::connect(eegwindow, &eegWindow::sendPrepStates, preProcessingworker, &preProcessingWorker::setPreprocessingParameters);
    QObject::connect(eegwindow, &eegWindow::setRemoveBCG, preProcessingworker, &preProcessingWorker::setRemoveBCG);
    QObject::connect(eegwindow, &eegWindow::set_processing_pause, preProcessingworker, &preProcessingWorker::set_processing_pause);
    QObject::connect(preProcessingworker, &preProcessingWorker::savePreprocessingOutput, &handler, &dataHandler::savePreprocessingOutput);
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
        MRIwin = new mriWindow(handler, this);
        connect(MRIwin, &TMSwindow::destroyed, this, &MainWindow::resetTMSwinPointer);
        
        // Add connection for MRI ROI updates
        if (MRIwin->getGlWidget()) {
            connect(MRIwin->getGlWidget(), &OpenGLMRIWidget::ROIListChanged, 
                    this, &MainWindow::updateMRIROIs);
        }
    }
    MRIwin->show();
    MRIwin->raise();
    MRIwin->activateWindow();
}

void MainWindow::resetMRIwinPointer() {
    MRIwin = nullptr;  // Reset the pointer after the window is destroyed
}

bool MainWindow::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::Close)
    {
        QDockWidget *dockWidget = qobject_cast<QDockWidget*>(obj);
        if (dockWidget)
        {
            // Remove from the QVector
            signalViewers.removeOne(dockWidget);

            // Optionally, perform any additional cleanup
            // For example, delete the dock widget
            dockWidget->deleteLater();
        }
    }

    // Pass the event on to the parent class
    return QMainWindow::eventFilter(obj, event);
}

void MainWindow::addSignalViewer()
{
    // Create a new dock widget
    QDockWidget *dockWidget = new QDockWidget(this);

    // Use the index in the QVector as the viewer number
    int viewerIndex = signalViewers.size();
    QString viewerTitle = QString("Signal Viewer %1").arg(viewerIndex + 1);
    dockWidget->setWindowTitle(viewerTitle);

    // Create the custom title bar with the sources, passing the dockWidget
    CustomTitleBar *titleBar = new CustomTitleBar(viewerTitle, signalSources, dockWidget, dockWidget);
    dockWidget->setTitleBarWidget(titleBar);

    // Create the content of the dock widget
    MainGlWidget* mainglWidget = new MainGlWidget(dockWidget);
    mainglWidget->setDataHandler(&handler);

    // Set size policies to allow expansion
    dockWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    mainglWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    // Set minimum height if desired
    mainglWidget->setMinimumHeight(50);

    // Set the widget
    dockWidget->setWidget(mainglWidget);

    // Set the allowed dock areas
    dockWidget->setAllowedAreas(Qt::AllDockWidgetAreas);

    // Add the dock widget to the main window
    if (signalViewers.isEmpty())
    {
        // First dock widget; add it to fill the space
        addDockWidget(Qt::TopDockWidgetArea, dockWidget);
    }
    else
    {
        // Not the first dock widget; split the previous one vertically
        QDockWidget *previousDockWidget = signalViewers.last();
        splitDockWidget(previousDockWidget, dockWidget, Qt::Vertical);
        
        // Calculate equal heights for all dock widgets
        int totalHeight = height();
        int equalHeight = totalHeight / (signalViewers.size() + 1);
        QList<QDockWidget*> allDocks = signalViewers;
        allDocks.append(dockWidget);
        QList<int> heights(allDocks.size(), equalHeight);
        
        // Resize all docks to equal heights
        resizeDocks(allDocks, heights, Qt::Vertical);
    }

    // Install the event filter on the dock widget
    dockWidget->installEventFilter(this);

    // Store the dock widget in the QVector
    signalViewers.append(dockWidget);

    // Connect signals
    connect(titleBar, &CustomTitleBar::signalSourceChanged, mainglWidget, &MainGlWidget::setSignalSource);
    connect(titleBar, &CustomTitleBar::triggerAVisibilityChanged, mainglWidget, &MainGlWidget::setTriggerAVisible);
    connect(titleBar, &CustomTitleBar::triggerBVisibilityChanged, mainglWidget, &MainGlWidget::setTriggerBVisible);
    connect(titleBar, &CustomTitleBar::triggerOutVisibilityChanged, mainglWidget, &MainGlWidget::setTriggerOutVisible);
    connect(titleBar, &CustomTitleBar::viewerTypeChanged, mainglWidget, &MainGlWidget::setViewerType);
    
    // Automatically set the signal source based on the viewer index
    mainglWidget->setSignalSource(viewerIndex);
}

void MainWindow::updateEEGChannels(const std::vector<std::string>& channelNames) {
    signalSources.clear();
    
    // Add EEG channel names
    for (const auto& name : channelNames) {
        signalSources.append(QString::fromStdString(name));
    }

    // Add existing ROI signals if available
    if (MRIwin && MRIwin->getGlWidget()) {
        const auto& ROInames = MRIwin->getGlWidget()->getROINames();
        for (const auto& name : ROInames) {
            signalSources.append(QString::fromStdString(name));
        }
    }

    updateSignalViewers();
}

void MainWindow::updateMRIROIs() {
    if (!MRIwin || !MRIwin->getGlWidget()) {
        return;
    }

    // Get current ROI names
    const auto& ROInames = MRIwin->getGlWidget()->getROINames();
    
    // Set ROI names in dataHandler
    handler.setROINames(ROInames);

    // Remove any existing ROI signals (keep only EEG channels)
    int numEEGChannels = handler.getChannelNames().size();
    while (signalSources.size() > numEEGChannels) {
        signalSources.removeLast();
    }
    
    // Add updated ROI signals
    for (const auto& name : ROInames) {
        signalSources.append(QString::fromStdString(name));
    }

    updateSignalViewers();
}

// Private helper method to update viewers and property
void MainWindow::updateSignalViewers() {
    // Update the property
    setProperty("signalSources", QVariant::fromValue(signalSources));
    
    // Update existing signal viewers
    for (auto* viewer : signalViewers) {
        if (auto* titleBar = qobject_cast<CustomTitleBar*>(viewer->titleBarWidget())) {
            titleBar->updateSignalSources(signalSources);
        }
    }
}