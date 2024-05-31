#ifndef PROCESSINGWINDOW_H
#define PROCESSINGWINDOW_H

#include <QMainWindow>
#include "processingglwidget.h"
#include "../dataHandler/dataHandler.h"

namespace Ui {
class ProcessingWindow;
}

class ProcessingWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit ProcessingWindow(dataHandler &handler, volatile std::sig_atomic_t &processingWorkerRunning, QWidget *parent = nullptr);
    ~ProcessingWindow();

public slots:
    void updateData();

private:
    Ui::ProcessingWindow *ui;

    dataHandler &handler;

    Eigen::MatrixXd processed_data;

    volatile std::sig_atomic_t &processingWorkerRunning;
};

#endif // PROCESSINGWINDOW_H
