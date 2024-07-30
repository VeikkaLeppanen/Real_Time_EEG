#include "processingglwidget.h"

ProcessingGlWidget::ProcessingGlWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    QTimer *timer = new QTimer(this);
    matrixCapasity_ = 30000;
    n_channels_ = 0;
    dataMatrix_ = Eigen::MatrixXd::Zero(n_channels_, matrixCapasity_);

    connect(timer, &QTimer::timeout, this, &ProcessingGlWidget::updateGraph);
    timer->start(16); // Update approximately every 16 ms (60 FPS)
}

void ProcessingGlWidget::initializeGL()
{
    initializeOpenGLFunctions();
    glClearColor(0.0, 0.0, 0.0, 1.0); // Set the clearing color to black
}

void ProcessingGlWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h); // Set the viewport to cover the entire widget
}

void ProcessingGlWidget::paintGL()
{
    if (n_channels_ == 0) return;  // Ensure there is data to draw
    
    // Initializing positional parameters
    int windowHeight = height();
    int enabled_channel_count = n_channels_;
    int rowHeight = windowHeight / std::max(1, enabled_channel_count);

    // min max values for y axis
    Eigen::VectorXd min_coeffs = Eigen::VectorXd::Zero(n_channels_);
    Eigen::VectorXd max_coeffs = Eigen::VectorXd::Zero(n_channels_);
    
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw graphs for enabled channels
    int graph_index = 0;
    for (int row = 0; row < dataMatrix_.rows(); row++) {
        // Set viewport for this row
        glViewport(0, graph_index * rowHeight, width(), rowHeight);

        // Set up the projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);  // Set the coordinate system to cover [-1,1] in both axes

        // Switch back to modelview matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Prepare data and draw the line strip
        Eigen::VectorXd dataVector = dataMatrix_.row(n_channels_ - 1 - row);
        
        double minVal, maxVal;
        if(useCustomScale_) {
            minVal = min_scale_;
            maxVal = max_scale_;
        } else {
            minVal = dataVector.minCoeff();
            maxVal = dataVector.maxCoeff();
        }
        min_coeffs(row) = minVal;
        max_coeffs(row) = maxVal;

        glColor3f(1.0, 1.0, 1.0); // Set the color to white for the lines
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < dataVector.size(); ++i) {
            float x = (float)i / (dataVector.size() - 1) * 2.0f - 1.0f;
            float y = 0;
            if (maxVal != minVal) {
                y = (dataVector(i) - minVal) / (maxVal - minVal) * 2.0f - 1.0f;
            }
            glVertex2f(x, y);
        }
        glEnd();

        graph_index++;
    }
    
    // Draw a vertical line at a specific index of the last graph
    int specific_index = 750; // Replace this with your specific index

    // Calculate the x-coordinate in NDC
    float x = (float)specific_index / (matrixCapasity_ - 1) * 2.0f - 1.0f;

    // Set viewport for the last graph
    glViewport(0, 0, width(), rowHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);  // Set the coordinate system to cover [-1,1] in both axes
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Draw the vertical line
    glColor3f(1.0, 0.0, 0.0); // Set the color to red for the vertical line
    glBegin(GL_LINES);
    glVertex2f(x, -1.0f);
    glVertex2f(x, 1.0f);
    glEnd();

    // QPainter for text overlays
    QPainter painter(this);
    painter.setPen(Qt::red);
    painter.setFont(QFont("Arial", 10)); // Set font here

    graph_index = 0;
    for (int row = 0; row < dataMatrix_.rows(); row++) {
        // Draw channel name for each channel
        int yPos_name = windowHeight - (rowHeight * graph_index + rowHeight / 2 + 10);  // Adjust vertical position
        QString name = "Undefined";  // Default name if no channel name is available
        if (row < channelNames_.size()) {
            name = channelNames_.at(n_channels_ - 1 - row);
        }
        painter.drawText(10, yPos_name, name);

        // Draw channel scales
        if (draw_channel_scales) {
            int yPos_min = windowHeight - (rowHeight * graph_index + 10);
            int yPos_max = windowHeight - (rowHeight * graph_index + rowHeight - 10);
            
            QString min_value = QString::number(min_coeffs(row));
            QString max_value = QString::number(max_coeffs(row));
            painter.drawText(width() - 50, yPos_min, min_value);
            painter.drawText(width() - 50, yPos_max, max_value);
        }


        graph_index++;
    }

    painter.end();
}

void ProcessingGlWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    emit fetchData();
    update();  // Request a re-draw
}
