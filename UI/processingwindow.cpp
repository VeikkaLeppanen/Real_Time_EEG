#include "processingwindow.h"
#include "ui_processingwindow.h"
#include <QThread>

ProcessingWindow::ProcessingWindow(dataHandler &handler, 
                    volatile std::sig_atomic_t &processingWorkerRunning, 
                               Eigen::MatrixXd &processed_data, 
                                       QWidget *parent)
    : QMainWindow(parent), 
      ui(new Ui::ProcessingWindow),
      handler(handler),
      processingWorkerRunning(processingWorkerRunning),
      processed_data(processed_data)
{
    ui->setupUi(this);

    // Initialize values for lineEdits
    ui->numberOfSamples->setText(QString::number(params.numberOfSamples));
    ui->downsampling->setText(QString::number(params.downsampling_factor));
    ui->delay->setText(QString::number(params.delay));
    ui->edge->setText(QString::number(params.edge));
    ui->modelOrder->setText(QString::number(params.modelOrder));
    ui->hilbertLength->setText(QString::number(params.hilbertWinLength));
    ui->stimulationTarget->setText(QString::number(params.stimulation_target));
    ui->phaseShift->setText(QString::number(params.phase_shift));

    setWindowTitle("Processing Window");
    resize(1280, 720);

    ProcessingGlWidget* processingglWidget = ui->processingGlWidget;
    if (processingglWidget) {

        connect(processingglWidget, &ProcessingGlWidget::fetchData, this, &ProcessingWindow::updateData);
        connect(this, &ProcessingWindow::setCustomScaleStatus, processingglWidget, &ProcessingGlWidget::setCustomScaleStatus);
        connect(this, &ProcessingWindow::setCustomScaleMin, processingglWidget, &ProcessingGlWidget::setCustomScaleMin);
        connect(this, &ProcessingWindow::setCustomScaleMax, processingglWidget, &ProcessingGlWidget::setCustomScaleMax);

        connect(this, &ProcessingWindow::getCustomScaleStatus, processingglWidget, &ProcessingGlWidget::getCustomScaleStatus);
        connect(this, &ProcessingWindow::getCustomScaleMin, processingglWidget, &ProcessingGlWidget::getCustomScaleMin);
        connect(this, &ProcessingWindow::getCustomScaleMax, processingglWidget, &ProcessingGlWidget::getCustomScaleMax);

        connect(this, &ProcessingWindow::updateWidgetChannelNames, processingglWidget, &ProcessingGlWidget::updateChannelNamesSTD);

        ui->checkBox->setChecked(emit getCustomScaleStatus());
        ui->scaleMin->setText(QString::number(emit getCustomScaleMin()));
        ui->scaleMax->setText(QString::number(emit getCustomScaleMax()));
    } else {
        // Error handling if glWidget is not found
        qWarning("Glwidget not found in UI!");
    }
}

void ProcessingWindow::updateData()
{
    ProcessingGlWidget* processingglWidget = ui->processingGlWidget;
    if (processingglWidget && processingWorkerRunning && (processed_data.size() > 0)) {

        std::lock_guard<std::mutex> lock(this->dataMutex); // Protect shared data access
        processingglWidget->updateMatrix(processed_data);
    }
}

ProcessingWindow::~ProcessingWindow()
{
    delete ui;
}

void ProcessingWindow::on_startButton_clicked()
{
    if(processingWorkerRunning) {
        QMessageBox::warning(this, "Error", "Processing is already running.");
    } else {
        std::cout << "ProcessingWindow start" << '\n';
        emit startProcessing(params);
    }
}


void ProcessingWindow::on_stopButton_clicked()
{
    processingWorkerRunning = 0;
}

void ProcessingWindow::on_numberOfSamples_editingFinished()
{
    bool ok;
    int value = ui->numberOfSamples->text().toInt(&ok);
    if (ok) {
        params.numberOfSamples = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}

void ProcessingWindow::on_downsampling_editingFinished()
{
    bool ok;
    int value = ui->downsampling->text().toInt(&ok);
    if (ok && (value > 0)) {
        params.downsampling_factor = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_delay_editingFinished()
{
    bool ok;
    int value = ui->delay->text().toInt(&ok);
    if (ok) {
        params.delay = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_edge_editingFinished()
{
    bool ok;
    int value = ui->edge->text().toInt(&ok);
    if (ok) {
        params.edge = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_modelOrder_editingFinished()
{
    bool ok;
    int value = ui->modelOrder->text().toInt(&ok);
    if (ok) {
        params.modelOrder = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_hilbertLength_editingFinished()
{
    bool ok;
    int value = ui->hilbertLength->text().toInt(&ok);
    if (ok) {
        params.hilbertWinLength = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_stimulationTarget_editingFinished()
{
    bool ok;
    double value = ui->stimulationTarget->text().toDouble(&ok);
    if (ok) {
        params.stimulation_target = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_phaseShift_editingFinished()
{
    bool ok;
    int value = ui->phaseShift->text().toInt(&ok);
    if (ok) {
        params.phase_shift = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}




void ProcessingWindow::on_scaleMax_editingFinished()
{
    bool ok;
    double value = ui->scaleMax->text().toDouble(&ok);
    if (ok) {
        emit setCustomScaleMax(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_scaleMin_editingFinished()
{
    bool ok;
    double value = ui->scaleMin->text().toDouble(&ok);
    if (ok) {
        emit setCustomScaleMin(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_checkBox_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setCustomScaleStatus(isChecked);
}

