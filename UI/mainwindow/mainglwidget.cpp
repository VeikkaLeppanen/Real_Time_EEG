#include "mainglwidget.h"

MainGlWidget::MainGlWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    windowLength_seconds = 2;   // Total time in seconds displayed
    time_line_spacing = 500;    // Line spacing in milliseconds
    QTimer *timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &MainGlWidget::updateGraph);
    timer->start(16); // Update approximately every 16 ms (60 FPS)
}

void MainGlWidget::initializeGL()
{
    initializeOpenGLFunctions();
    glClearColor(0.0, 0.0, 0.0, 1.0); // Set the clearing color to black
}

void MainGlWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h); // Set the viewport to cover the entire widget
}

void MainGlWidget::paintGL()
{

}


void MainGlWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    update();  // Request a re-draw
}

void MainGlWidget::setSignalSource(const QString &source)
{
    // Update the signal source based on the selection
    // Implement the logic to handle the signal source change
    qDebug() << "Signal source changed to:" << source;

    // For example, you might have a method like:
    // loadSignalSource(source);
}

void MainGlWidget::updateData(const Eigen::VectorXi &newData, 
                              const Eigen::VectorXi &triggers_A, 
                              const Eigen::VectorXi &triggers_B, 
                              const Eigen::VectorXd &time_stamps, 
                                        std::string source_name) { 
    if (!pause_view) {
        std::lock_guard<std::mutex> lock(this->dataMutex);
        
        // Check if dimensions differ
        if (dataVector_.size() != newData.size()) dataVector_.resize(newData.size());
        if (triggers_A_.size() != triggers_A.size()) triggers_A_.resize(triggers_A.size());
        if (triggers_B_.size() != triggers_B.size()) triggers_B_.resize(triggers_B.size());
        if (time_stamps_.size() != time_stamps.size()) time_stamps_.resize(time_stamps.size());

        dataVector_ = newData;
        vectorCapasity = newData.size();
        triggers_A_ = triggers_A;
        triggers_B_ = triggers_B;
        time_stamps_ = time_stamps;
        source_name_ = QString::fromStdString(source_name);
    }
}