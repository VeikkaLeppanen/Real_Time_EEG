#include "processingwindow.h"
#include "ui_processingwindow.h"

ProcessingWindow::ProcessingWindow(dataHandler &handler, volatile std::sig_atomic_t &processingWorkerRunning, QWidget *parent)
    : QMainWindow(parent), 
      ui(new Ui::ProcessingWindow),
      handler(handler),
      processingWorkerRunning(processingWorkerRunning)
{
    ui->setupUi(this);

    setWindowTitle("Processing Window");

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
    if (processingglWidget/* && processingWorkerRunning && (processed_data.size() > 0)*/) {

        processingglWidget->updateMatrix(processed_data);
    }
}

ProcessingWindow::~ProcessingWindow()
{
    delete ui;
}
