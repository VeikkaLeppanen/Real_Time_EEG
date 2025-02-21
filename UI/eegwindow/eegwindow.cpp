#include "eegwindow.h"

eegWindow::eegWindow(dataHandler &handler, 
      volatile std::sig_atomic_t &signal_received,
      volatile std::sig_atomic_t &processingWorkerRunning,
                         QWidget *parent) 
    : QMainWindow(parent), 
      ui(new Ui::EegWindow),
      handler(handler),
      signal_received(signal_received),
      processingWorkerRunning(processingWorkerRunning)
    {
        try {
            ui->setupUi(this);

            // Setting validators for lineEdits
            ui->lineEditPort->setValidator(new QIntValidator(0, 65535, this));
            ui->lineEditTimeOut->setValidator(new QIntValidator(0, 1000, this));
            ui->lineEditTALength->setValidator(new QIntValidator(0, 500000, this));
            ui->lineEditTRLength->setValidator(new QIntValidator(0, 500000, this));
            ui->lineEditGAaverage->setValidator(new QIntValidator(0, 1000, this));

            // Initialize values for lineEdits
            ui->lineEditPort->setText(QString::number(port));  // Example port number
            ui->lineEditTimeOut->setText(QString::number(bridge_timeout)); // Example timeout in milliseconds
            ui->lineEditTALength->setText(QString::number(handler.get_TA_length()));
            ui->lineEditTRLength->setText(QString::number(handler.get_TR_length()));
            ui->lineEditGAaverage->setText(QString::number(handler.get_GA_average_length()));
            ui->lineEdit_XaxisSpacing->setText(QString::number(500));
            ui->checkBox_GA->setChecked(handler.getGAState());
            ui->checkBox_2->setChecked(handler.getBaselineState());
            ui->filter1->setChecked(handler.getFilterState());
            ui->numberOfSamples->setText(QString::number(prepParams.numberOfSamples));
            ui->downsampling->setText(QString::number(prepParams.downsampling_factor));
            ui->delay->setText(QString::number(prepParams.delay));

            ui->tabWidget->setTabText(0, "Device");
            ui->tabWidget->setTabText(1, "Preprocessing");
            ui->checkBox_triggers_A->setStyleSheet("QCheckBox { color : blue; }");
            ui->checkBox_triggers_B->setStyleSheet("QCheckBox { color : green; }");
            setWindowTitle("EEG Window");
            resize(1280, 720);

            checkHandlerTimer = new QTimer(this);
            connect(checkHandlerTimer, &QTimer::timeout, this, &eegWindow::checkHandlerReady);

            glWidget = ui->openglWidget;
            if (glWidget) {
                connect(glWidget, &Glwidget::fetchData, this, &eegWindow::updateData);
                connect(this, &eegWindow::updateChannelNamesSTD, glWidget, &Glwidget::updateChannelNamesSTD);
                connect(this, &eegWindow::updateChannelNamesQt, glWidget, &Glwidget::updateChannelNamesQt);
                connect(this, &eegWindow::updateChannelDisplayState, glWidget, &Glwidget::updateChannelDisplayState);
                connect(this, &eegWindow::scaleDrawStateChanged, glWidget, &Glwidget::scaleDrawStateChanged);
                connect(this, &eegWindow::setShowTriggers_A, glWidget, &Glwidget::setShowTriggers_A);
                connect(this, &eegWindow::setShowTriggers_B, glWidget, &Glwidget::setShowTriggers_B);
                connect(this, &eegWindow::switchPause, glWidget, &Glwidget::switchPause);
                connect(this, &eegWindow::setDrawXaxis, glWidget, &Glwidget::setDrawXaxis);
                connect(this, &eegWindow::updateTLineSpacing, glWidget, &Glwidget::updateTLineSpacing);

                ui->checkBox_triggers_A->setChecked(glWidget->getShowTriggers_A());
                ui->checkBox_triggers_B->setChecked(glWidget->getShowTriggers_B());
                ui->checkBox_Xaxis->setChecked(glWidget->getDrawXaxis());

                emit updateChannelNamesSTD(handler.getChannelNames());
            } else {
                // Error handling if glWidget is not found
                qWarning("Glwidget not found in UI!");
            }
        } catch (std::exception& e) {
            emit error(QString("An error occurred during eegWindow initialization: %1").arg(e.what()));
        }
    }

eegWindow::~eegWindow() {
    delete ui;
}

void eegWindow::newBCGState(bool state) {
    ui->checkBox_removeBCG->setChecked(state);
}

void eegWindow::handleError(const QString &error)
{
    // Error handling code here
    QMessageBox::critical(this, "Error", error);
}

void eegWindow::updateData() 
{
    if (processingWorkerRunning) return;

    glWidget = ui->openglWidget;
    if (glWidget && handler.isReady()) {
        Eigen::MatrixXd data;
        Eigen::VectorXi triggers_A;
        Eigen::VectorXi triggers_B;
        Eigen::VectorXi triggers_out;
        Eigen::VectorXd time_stamps;
        handler.getLatestDataAndTriggers(data, triggers_A, triggers_B, triggers_out, time_stamps, samples_to_display);
        glWidget->updateMatrix(data, triggers_A, triggers_B, time_stamps, handler.getChannelNames());
    }
}

void eegWindow::on_connectButton_clicked()
{
    emit connectEegBridge(port, bridge_timeout);

    checkHandlerTimer->start(100);  // Check every 500 ms
}

void eegWindow::checkHandlerReady() {
    if (handler.isReady()) {
        checkHandlerTimer->stop();  // Stop the timer as we don't need to check anymore

        source_channels_ = handler.getSourceChannels();

        if (channelMap_.size() > 0) {
            setupChannelNames();
        } else {
            QStringList QchannelNames;
            for (size_t i = 0; i < source_channels_.size(); i++) {
                QchannelNames.append("Undefined");
            }

            channelNames_ = QchannelNames;

            emit updateChannelNamesQt(QchannelNames);
        }
        preProcessing_start();
        emit requestBCGState();
    }
}


void eegWindow::on_disconnectButton_clicked()
{
    handler.setHandlerState(WAITING_FOR_START);
    preProcessing_stop();
    signal_received = 1;
}

void eegWindow::on_sourceChannelLoad_clicked()
{
    std::map<int, std::string> channelMapStd;
    QString fileName = QFileDialog::getOpenFileName(
        this,                 // parent widget
        "Open Document",      // dialog caption
        QDir::homePath(),     // starting directory
        "Text Files (*.txt)", // file types
        nullptr,
        QFileDialog::DontUseNativeDialog
    );

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QTextStream in(&file);
            QString line = in.readLine();  // Read and discard the header line, if any

            while (!in.atEnd()) {
                line = in.readLine().trimmed();
                if (line.isEmpty()) {
                    continue; // Skip empty lines
                }

                QStringList fields = line.split('\t');  // assuming tab-separated values
                if (fields.size() >= 2) {  // Check if there are at least two fields: key and value
                    bool ok;
                    int key = fields[0].toInt(&ok);
                    if (ok) {
                        std::string value = fields[1].toStdString();
                        channelMapStd[key] = value;
                    } else {
                        qDebug() << "Invalid key in line: " << line;
                    }
                } else {
                    qDebug() << "Invalid line format: " << line;
                }
            }
            file.close();
        } else {
            // Handle error if the file can't be opened
            qDebug("Failed to open the file for reading.");
        }

        channelMap_ = channelMapStd;

        if (source_channels_.size() > 0) {
            setupChannelNames();
        }
        std::cout << "Channel map loaded successfully." << std::endl;
    }
}

void eegWindow::setupChannelNames()
{
    std::vector<std::string> channelNames;
    QStringList QchannelNames;

    for (size_t i = 0; i < source_channels_.size(); i++) {
        int channelIndex = source_channels_(i);

        // Check if the channel index exists in the map
        auto it = channelMap_.find(channelIndex);
        if (it != channelMap_.end()) {
            // Key exists, safe to use
            std::string C_Name = it->second;
            channelNames.push_back(C_Name);
            QchannelNames.append(QString::fromStdString(C_Name));
        } else {
            // Handle the case where the channel index is not found
            qDebug() << "Channel index" << channelIndex << "not found in channelMap_.";

            // Option 1: Assign a default name
            std::string defaultName = "Undefined";
            channelNames.push_back(defaultName);
            QchannelNames.append(QString::fromStdString(defaultName));
        }
    }

    handler.setChannelNames(channelNames);

    channelNames_ = QchannelNames;

    emit updateChannelNamesQt(QchannelNames);
}

void eegWindow::on_lineEditPort_editingFinished()
{
    bool ok;
    int value = ui->lineEditPort->text().toInt(&ok);
    if (ok) {
        port = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 65535.");
    }
}

void eegWindow::on_lineEditTimeOut_editingFinished()
{
    bool ok;
    int value = ui->lineEditTimeOut->text().toInt(&ok);
    if (ok) {
        bridge_timeout = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 1000.");
    }
}

void eegWindow::on_lineEditTALength_editingFinished()
{
    bool ok1, ok2;
    int value1 = ui->lineEditTALength->text().toInt(&ok1);
    int value2 = ui->lineEditGAaverage->text().toInt(&ok2);
    if (ok1 && ok2) {
        emit stopGACorrection();
        emit applyGACorrection(value1, value2);
        emit startGACorrection();
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 500000.");
    }
}

void eegWindow::on_lineEditTRLength_editingFinished()
{
    bool ok1, ok2, ok3;
    int value1 = ui->lineEditTRLength->text().toInt(&ok1);
    int value2 = ui->lineEditTALength->text().toInt(&ok2);
    int value3 = ui->lineEditGAaverage->text().toInt(&ok3);
    if (ok1 && ok2 && ok3) {
        handler.set_TR_length(value1);
        emit stopGACorrection();
        emit applyGACorrection(value2, value3);
        emit startGACorrection();
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 500000.");
    }
}

void eegWindow::on_lineEditGAaverage_editingFinished()
{
    bool ok1, ok2;
    int value1 = ui->lineEditTALength->text().toInt(&ok1);
    int value2 = ui->lineEditGAaverage->text().toInt(&ok2);
    if (ok1 && ok2) {
        emit stopGACorrection();
        emit applyGACorrection(value1, value2);
        emit startGACorrection();
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 1000.");
    }
}


void eegWindow::on_checkBox_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit scaleDrawStateChanged(isChecked);
}


void eegWindow::on_filter1_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    handler.setFilterState(isChecked);
}


void eegWindow::on_checkBox_2_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    handler.setBaselineState(isChecked);
}


void eegWindow::on_numberOfSamples_editingFinished()
{
    bool ok;
    int value = ui->numberOfSamples->text().toInt(&ok);
    if (ok) {
        if (prepParams.numberOfSamples != value) {
            prepParams.numberOfSamples = value;
            updateChannelLength(value);
            
            emit set_processing_pause(true);
            emit sendPrepStates(prepParams);
            emit set_processing_pause(false);
        }
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}

void eegWindow::on_downsampling_editingFinished()
{
    bool ok;
    int value = ui->downsampling->text().toInt(&ok);
    if (ok) {
        if (prepParams.downsampling_factor != value) {
            prepParams.downsampling_factor = value;
            
            emit set_processing_pause(true);
            emit sendPrepStates(prepParams);
            emit set_processing_pause(false);
        }
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void eegWindow::on_delay_editingFinished()
{
    bool ok;
    int value = ui->delay->text().toInt(&ok);
    if (ok) {
        if (prepParams.delay != value) {
            prepParams.delay = value;

            emit set_processing_pause(true);
            emit sendPrepStates(prepParams);
            emit set_processing_pause(false);
        }
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}

void eegWindow::preProcessing_start()
{
    if(processingWorkerRunning) {
        QMessageBox::warning(this, "Error", "Processing is already running.");
    } else {
        std::cout << "Preprocessing start" << '\n';
        emit startPreprocessing(prepParams, phaseEstParams);
        preprocessingWorkerRunning = true;

        bool ok;
        int value = ui->numberOfSamples->text().toInt(&ok);
        if (ok) {
            samples_to_display = value;
            updateChannelLength(value);
        } else {
            QMessageBox::warning(this, "Window Length Error", " ");
        }
    }
}

void eegWindow::preProcessing_stop()
{
    std::cout << "Preprocessing stop" << '\n';
    preprocessingWorkerRunning = false;
    processingWorkerRunning = 0;
}

void eegWindow::on_checkBox_GA_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    bool ok1, ok2;
    int value1 = ui->lineEditTALength->text().toInt(&ok1);
    int value2 = ui->lineEditGAaverage->text().toInt(&ok2);
    if (isChecked && ok1 && ok2) {
        emit stopGACorrection();
        emit applyGACorrection(value1, value2);
        emit startGACorrection();
    } else {
        emit stopGACorrection();
    }
}


void eegWindow::on_checkBox_triggers_A_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_A(isChecked);
}


void eegWindow::on_checkBox_triggers_B_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_B(isChecked);
}


void eegWindow::on_pushButton_pauseView_clicked()
{
    emit switchPause();
}


void eegWindow::on_checkBox_removeBCG_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setRemoveBCG(isChecked);
    emit set_processing_pause(true);
    emit sendPrepStates(prepParams);
    emit set_processing_pause(false);
}

void eegWindow::on_checkBox_Xaxis_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setDrawXaxis(isChecked);
}

void eegWindow::on_lineEdit_XaxisSpacing_editingFinished()
{
    bool ok;
    int value = ui->lineEdit_XaxisSpacing->text().toInt(&ok);
    if (ok) {
        emit updateTLineSpacing(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}

void eegWindow::updateChannelLength(int value)
{
    if (glWidget) {
        glWidget->updateWindowLength_seconds(value / handler.getSamplingRate());
    }
}
void eegWindow::on_checkBox_GeoSum_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    handler.setGeometricSumState(isChecked);
}

