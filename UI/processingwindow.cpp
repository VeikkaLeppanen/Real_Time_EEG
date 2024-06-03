#include "processingwindow.h"
#include "ui_processingwindow.h"
#include <QThread>

ProcessingWindow::ProcessingWindow(dataHandler &handler, volatile std::sig_atomic_t &processingWorkerRunning, Eigen::MatrixXd &processed_data, QWidget *parent)
    : QMainWindow(parent), 
      ui(new Ui::ProcessingWindow),
      handler(handler),
      processingWorkerRunning(processingWorkerRunning),
      processed_data(processed_data)
{
    ui->setupUi(this);

    // Initialize values for lineEdits
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

    } else {
        // Error handling if glWidget is not found
        qWarning("Glwidget not found in UI!");
    }
}

void ProcessingWindow::updateData()
{
    ProcessingGlWidget* processingglWidget = ui->processingGlWidget;
    if (processingglWidget && processingWorkerRunning && (processed_data.size() > 0)) {

        processingglWidget->updateMatrix(processed_data);
    }
}

ProcessingWindow::~ProcessingWindow()
{
    delete ui;
}

void ProcessingWindow::on_startButton_clicked()
{
    processingWorkerRunning = 0;
    QThread::msleep(10);
    std::cout << "ProcessingWindow start" << '\n';
    emit startProcessing(params);
}


void ProcessingWindow::on_stopButton_clicked()
{
    processingWorkerRunning = 0;
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

