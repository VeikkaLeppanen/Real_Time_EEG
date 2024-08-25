#include "processingglwidget.h"

ProcessingGlWidget::ProcessingGlWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    QTimer *timer = new QTimer(this);
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
    if (dataMatrix_.rows() == 0) return;
    
    // Initializing positional parameters
    int windowHeight = height();
    int rowHeight = windowHeight / n_channels_;

    // min max values for y axis
    Eigen::VectorXd min_coeffs = Eigen::VectorXd::Zero(n_channels_);
    Eigen::VectorXd max_coeffs = Eigen::VectorXd::Zero(n_channels_);
    
    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_LINE_STIPPLE);

    float currentTimePosition = static_cast<float>(numPastElements_) / (numPastElements_ + numFutureElements_);
    float glX = currentTimePosition * 2.0f - 1.0f; // Convert to OpenGL's coordinate system

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
    
    // Set viewport to cover the entire widget for drawing triggers
    glViewport(0, 0, width(), windowHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    
    // Draw triggers over all graphs
    int totalDataPoints = triggers_A_.size();
    if (show_triggers_A) {
        for (int i = 0; i < triggers_A_.size(); ++i) {
            if (triggers_A_(i) == 1) {
                float x = ((float)i / (totalDataPoints - 1)) * (glX + 1) - 1;
                glColor3f(0.0, 0.0, 1.0);
                glBegin(GL_LINES);
                glVertex2f(x, -1.0);
                glVertex2f(x, 1.0);
                glEnd();
            }
        }
    }

    if (show_triggers_B) {
        for (int i = 0; i < triggers_B_.size(); ++i) {
            if (triggers_B_(i) == 1) {
                float x = ((float)i / (totalDataPoints - 1)) * (glX + 1) - 1;
                glColor3f(0.0, 1.0, 0.0);
                glBegin(GL_LINES);
                glVertex2f(x, -1.0);
                glVertex2f(x, 1.0);
                glEnd();
            }
        }
    }
    
    // Enable and configure line stipple for dashed line
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x00FF);  // 1x repeat factor, 0x00FF pattern

    // Draw the current time line across the full height
    glColor3f(1.0, 0.0, 0.0); // Red for the current time line
    glBegin(GL_LINES);
    glVertex2f(glX, -1.0);
    glVertex2f(glX, 1.0);
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

void ProcessingGlWidget::updateMatrix(const Eigen::MatrixXd &newMatrix, 
                                      const Eigen::VectorXi &triggers_A, 
                                      const Eigen::VectorXi &triggers_B, 
                                                        int numPastElements, 
                                                        int numFutureElements) 
{
    if (!pause_view) {
        // Check if dimensions differ
        if (dataMatrix_.rows() != newMatrix.rows() || dataMatrix_.cols() != newMatrix.cols()) {
            // Resize dataMatrix_ to match the dimensions of newMatrix
            dataMatrix_.resize(newMatrix.rows(), newMatrix.cols());
        }

        // Assign the new matrix
        dataMatrix_ = newMatrix;

        // Update other attributes based on the new matrix
        matrixCapasity_ = newMatrix.cols(); 
        n_channels_ = newMatrix.rows();
        numPastElements_ = numPastElements;
        numFutureElements_ = numFutureElements;
        triggers_A_ = triggers_A;
        triggers_B_ = triggers_B;
    }
}

void ProcessingGlWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    // emit fetchData();
    update();  // Request a re-draw
}

void ProcessingGlWidget::updateChannelNamesSTD(std::vector<std::string> channelNames) {
    QStringList newNames;
    for(size_t i = 0; i < channelNames.size(); i++) {
        newNames.append(QString::fromStdString(channelNames[i])); 
    }
    channelNames_ = newNames;
}