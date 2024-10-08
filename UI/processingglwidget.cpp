#include "processingglwidget.h"

ProcessingGlWidget::ProcessingGlWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    windowLength_seconds = 2;   // Total time in seconds displayed
    time_line_spacing = 500;    // Line spacing in milliseconds
    totalTimeLines = (windowLength_seconds * 1000) / time_line_spacing;  // Calculate number of lines to draw
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
    int bottomMarging = 20;
    int windowHeight = height() - bottomMarging;
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
        glViewport(0, bottomMarging + graph_index * rowHeight, width(), rowHeight);

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
    glViewport(0, 0, width(), windowHeight + bottomMarging);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    
    // Draw trigger lines for triggers_A
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

    // Draw trigger lines for triggers_B
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

    // Draw trigger lines for triggers_out
    if (show_triggers_out) {
        for (int i = 0; i < triggers_out_.size(); ++i) {
            if (triggers_out_(i) == 1) {
                float x = ((float)i / (totalDataPoints - 1)) * (glX + 1) - 1;
                glColor3f(0.0, 1.0, 1.0);
                glBegin(GL_LINES);
                glVertex2f(x, -1.0);
                glVertex2f(x, 1.0);
                glEnd();
            }
        }
    }
    
    // Draw the vertical line for the current time
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x00FF);  // 1x repeat factor, 0x00FF pattern
    glColor3f(1.0, 0.0, 0.0); // Red for the current time line
    glBegin(GL_LINES);
    glVertex2f(glX, -1.0);
    glVertex2f(glX, 1.0);
    glEnd();

    // Initialize tracker to the next whole second
    double time_line_spacing_seconds = time_line_spacing * 1e-3;
    double initialTimeInSeconds = time_stamps_(0) / 1e3;
    double tracker = std::ceil(initialTimeInSeconds / time_line_spacing_seconds) * time_line_spacing_seconds;

    if (drawXaxis) {
        // Enable stipple for dashed lines
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1, 0x00FF);  // 1x repeat factor, 0x00FF pattern
        glColor3f(1.0, 0.0, 0.0);  // Gray color for the lines

        // Draw each vertical line
        for (int i = 0; i < totalDataPoints; i++) {
            double timeInSeconds = time_stamps_(i) / 1e3; // Convert microseconds to seconds
            if (timeInSeconds >= tracker) {
                float x = ((float)i / (totalDataPoints - 1)) * (glX + 1) - 1;  // Convert index to OpenGL coordinates
                glBegin(GL_LINES);
                glVertex2f(x, -1.0);
                glVertex2f(x, 1.0);
                glEnd();
                tracker += time_line_spacing_seconds;
            }
        }

        glDisable(GL_LINE_STIPPLE); // Disable stipple after drawing lines
    }

    // Draw the channel names
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

    tracker = std::ceil(initialTimeInSeconds);  // Round up to next whole second
    painter.setFont(QFont("Arial", 18)); // Set font here

    // Draw the timestamp labels and vertical lines
    for (int i = 0; i < totalDataPoints; ++i) {
        double timeInSeconds = time_stamps_(i) / 1e3; // Convert microseconds to seconds
        if (timeInSeconds >= tracker) {
            QString label = QString::number(tracker, 'f', 0); // Label for whole seconds
            float x = (float)i / totalDataPoints * (width() * currentTimePosition); // Calculate pixel x-coordinate

            // Draw text label
            painter.drawText(x, windowHeight + bottomMarging, label); // Draw just above the bottom

            tracker += 1.0;  // Increment tracker by one second
        }
    }

    painter.end();
}

void ProcessingGlWidget::updateMatrix(const Eigen::MatrixXd &newMatrix, 
                                      const Eigen::VectorXi &triggers_A, 
                                      const Eigen::VectorXi &triggers_B, 
                                      const Eigen::VectorXi &triggers_out, 
                                      const Eigen::VectorXd &time_stamps,
                                                        int numPastElements, 
                                                        int numFutureElements) 
{
    if (!pause_view) {
        std::lock_guard<std::mutex> lock(this->dataMutex);

        // Check if dimensions differ
        if (dataMatrix_.rows() != newMatrix.rows() || dataMatrix_.cols() != newMatrix.cols()) dataMatrix_.resize(newMatrix.rows(), newMatrix.cols());
        if (triggers_A_.size() != triggers_A.size()) triggers_A_.resize(triggers_A.size());
        if (triggers_B_.size() != triggers_B.size()) triggers_B_.resize(triggers_B.size());
        if (triggers_out_.size() != triggers_out.size()) triggers_out_.resize(triggers_out.size());
        if (time_stamps_.size() != time_stamps.size()) time_stamps_.resize(time_stamps.size());

        // Assign the new matrix
        dataMatrix_ = newMatrix;

        // Update other attributes based on the new matrix
        matrixCapasity_ = newMatrix.cols(); 
        n_channels_ = newMatrix.rows();
        numPastElements_ = numPastElements;
        numFutureElements_ = numFutureElements;
        triggers_A_ = triggers_A;
        triggers_B_ = triggers_B;
        triggers_out_ = triggers_out;
        time_stamps_ = time_stamps;
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