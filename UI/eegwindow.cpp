#include "eegwindow.h"

#include <QLabel>

eegWindow::eegWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent) 
    : QMainWindow(parent), 
      ui(new Ui::EegWindow),
      handler(handler),
      signal_received(signal_received)
    {
        ui->setupUi(this);
        // QLabel *label = new QLabel("Hello from the second window!", this);
        setWindowTitle("EEG Window");
        // Optional: resize(800, 600); // Set a more appropriate default size

        Glwidget* glWidget = ui->openglWidget;  // Assuming the object name in the .ui is glWidget
        if (glWidget) {
            connect(glWidget, &Glwidget::fetchData, this, &eegWindow::updateData);
        } else {
            // Error handling if glWidget is not found
            qWarning("Glwidget not found in UI!");
        }
    }

eegWindow::~eegWindow() {
    delete ui;  // Clean up the ui pointer
}

void eegWindow::updateData() 
{
    Glwidget* glWidget = ui->openglWidget;  // Assuming ui->glWidget is the pointer to the Glwidget
    if (glWidget && handler.isReady()) {
        Eigen::MatrixXd newMatrix = handler.returnLatestDataInOrder(samples_to_display);
    
        glWidget->updateMatrix(newMatrix);  // You might need to implement updateMatrix in Glwidget
    }
}

void eegWindow::on_connectButton_clicked()
{
    if (!bridge.isRunning()) {
        signal_received = 0;

        bridge.setPort(port);
        handler.seg_GAcorr_params(GALength, GAAverage);

        QThread* thread = new QThread;
        Worker* worker = new Worker(bridge, handler, signal_received);
        worker->moveToThread(thread);

        connect(thread, &QThread::started, worker, &Worker::process);
        connect(worker, &Worker::finished, thread, &QThread::quit);
        connect(worker, &Worker::error, this, &eegWindow::handleError); // Handle errors
        connect(worker, &Worker::finished, worker, &Worker::deleteLater);
        connect(thread, &QThread::finished, thread, &QThread::deleteLater);

        thread->start();
    }
}

void eegWindow::handleError(const QString &error)
{
    // Error handling code here
    QMessageBox::critical(this, "Error", error);
}

void eegWindow::on_disconnectButton_clicked()
{
    signal_received = 1;
}

