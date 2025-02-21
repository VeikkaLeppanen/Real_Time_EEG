#include "mainglwidget.h"

MainGlWidget::MainGlWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    windowLength_seconds = 2;   // Total time in seconds displayed
    time_line_spacing = 500;    // Line spacing in milliseconds
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
    // Don't draw if there's no data
    if (dataVector_.size() == 0) return;

    // Call appropriate painting function based on viewer type
    if (viewerType == 0) {
        glClear(GL_COLOR_BUFFER_BIT);
        paintSignalViewer();
    } else {
        bool hasValidTrigger = false;
        
        // Calculate trigger position with wraparound using modulo
        int bufferSize = dataVector_.size();
        int triggerPosition = ((data_index_ - response_window_size / 2) + bufferSize) % bufferSize;
        
        // Check for triggers within the tolerance range
        for (int i = -trigger_position_tolerance; i <= trigger_position_tolerance; i++) {
            int checkPos = ((triggerPosition + i) + bufferSize) % bufferSize;
            if (triggers_out_[checkPos] != 0) {
                hasValidTrigger = true;
                break;
            }
        }
        
        if (hasValidTrigger) {
            glClear(GL_COLOR_BUFFER_BIT);
            paintResponseMonitor();
        }
    }
}

void MainGlWidget::paintSignalViewer()
{
    // Lock data while drawing
    std::lock_guard<std::mutex> lock(this->dataMutex);

    // Set up the projection matrix for normalized coordinates
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Calculate min/max values
    auto [minVal, maxVal] = calculateDataRange(dataVector_);

    // Draw the signal
    glColor3f(1.0f, 1.0f, 1.0f);  // White for the signal
    glBegin(GL_LINE_STRIP);
    
    if (currentSignalType == SignalType::MRI) {
        // For MRI signals, draw step function with minimal vertices
        double currentValue = dataVector_[0];
        float y = 0;
        if (maxVal != minVal) {
            y = (currentValue - minVal) / (maxVal - minVal) * 2.0f - 1.0f;
        }
        
        // Add first point
        glVertex2f(-1.0f, y);
        
        for (int i = 1; i < dataVector_.size(); i++) {
            float x = (static_cast<float>(i) / (dataVector_.size() - 1)) * 2.0f - 1.0f;
            
            // If value changes, add two vertices to create the step
            if (dataVector_[i] != currentValue) {
                // Add vertex at previous y value (creates horizontal line)
                // glVertex2f(x, y);
                
                // Update to new value
                currentValue = dataVector_[i];
                if (maxVal != minVal) {
                    y = (currentValue - minVal) / (maxVal - minVal) * 2.0f - 1.0f;
                }
                
                // Add vertex at new y value (creates vertical line)
                glVertex2f(x, y);
            }
        }
        
        // Add final point at the right edge
        glVertex2f(1.0f, y);
    } else {
        // Original drawing code for EEG signals
        for (int i = 0; i < dataVector_.size(); i++) {
            float x = (static_cast<float>(i) / (dataVector_.size() - 1)) * 2.0f - 1.0f;
            float y = 0;
            if (maxVal != minVal) {
                y = (dataVector_[i] - minVal) / (maxVal - minVal) * 2.0f - 1.0f;
            }
            glVertex2f(x, y);
        }
    }
    glEnd();

    // Draw triggers as vertical lines
    // Trigger A (Blue)
    if (show_triggers_A) {
        glColor3f(0.0f, 0.0f, 1.0f);  // Blue
        glBegin(GL_LINES);
        for (int i = 0; i < triggers_A_.size(); i++) {
            if (triggers_A_[i] != 0) {
                float x = (static_cast<float>(i) / (dataVector_.size() - 1)) * 2.0f - 1.0f;
                glVertex2f(x, -1.0f);  // Bottom of screen
                glVertex2f(x, 1.0f);   // Top of screen
            }
        }
        glEnd();
    }

    // Trigger B (Green)
    if (show_triggers_B) {
        glColor3f(0.0f, 1.0f, 0.0f);  // Green
        glBegin(GL_LINES);
        for (int i = 0; i < triggers_B_.size(); i++) {
            if (triggers_B_[i] != 0) {
                float x = (static_cast<float>(i) / (dataVector_.size() - 1)) * 2.0f - 1.0f;
                glVertex2f(x, -1.0f);
                glVertex2f(x, 1.0f);
            }
        }
        glEnd();
    }

    // Trigger Out (Sky Blue) - assuming trigger out is when both A and B are active
    if (show_trigger_out) {
        glColor3f(0.529f, 0.808f, 0.922f);  // Sky Blue
        glBegin(GL_LINES);
        for (int i = 0; i < triggers_out_.size(); i++) {
            if (triggers_out_[i] != 0) {
                float x = (static_cast<float>(i) / (dataVector_.size() - 1)) * 2.0f - 1.0f;
                glVertex2f(x, -1.0f);
                glVertex2f(x, 1.0f);
            }
        }
        glEnd();
    }

    // Render min/max values using QPainter
    QPainter painter(this);
    painter.setPen(Qt::red);
    QFont font = painter.font();
    font.setPointSize(10);
    painter.setFont(font);

    QString maxText = QString::number(maxVal, 'f', 2);
    QString minText = QString::number(minVal, 'f', 2);

    painter.drawText(width() - 70, 20, maxText);
    painter.drawText(width() - 70, height() - 10, minText);
    
    painter.end();
}

void MainGlWidget::paintResponseMonitor()
{
    // Lock data while drawing
    std::lock_guard<std::mutex> lock(this->dataMutex);

    // Set up the projection matrix for normalized coordinates
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Calculate the start index with wraparound
    int bufferSize = dataVector_.size();
    int startIndex = ((data_index_ - response_window_size) + bufferSize) % bufferSize;

    // Create a temporary vector for the window of data
    Eigen::VectorXd windowData(response_window_size);
    for (int i = 0; i < response_window_size; i++) {
        int idx = (startIndex + i) % bufferSize;
        windowData[i] = dataVector_[idx];
    }

    // Calculate min/max values for the window
    auto [minVal, maxVal] = calculateDataRange(windowData);

    // Draw the signal
    glColor3f(1.0f, 1.0f, 1.0f);  // White for the signal
    glBegin(GL_LINE_STRIP);
    
    // Draw each sample in the window
    for (int i = 0; i < response_window_size; i++) {
        float x = (static_cast<float>(i) / (response_window_size - 1)) * 2.0f - 1.0f;
        
        float y = 0;
        if (maxVal != minVal) {
            y = (windowData[i] - minVal) / (maxVal - minVal) * 2.0f - 1.0f;
        }
        
        glVertex2f(x, y);
    }
    glEnd();

    // Draw triggers within the window
    for (int i = 0; i < response_window_size; i++) {
        int idx = (startIndex + i) % bufferSize;
        float x = (static_cast<float>(i) / (response_window_size - 1)) * 2.0f - 1.0f;

        // Trigger A (Blue)
        if (show_triggers_A && triggers_A_[idx] != 0) {
            glColor3f(0.0f, 0.0f, 1.0f);
            glBegin(GL_LINES);
            glVertex2f(x, -1.0f);
            glVertex2f(x, 1.0f);
            glEnd();
        }

        // Trigger B (Green)
        if (show_triggers_B && triggers_B_[idx] != 0) {
            glColor3f(0.0f, 1.0f, 0.0f);
            glBegin(GL_LINES);
            glVertex2f(x, -1.0f);
            glVertex2f(x, 1.0f);
            glEnd();
        }

        // Trigger Out (Sky Blue)
        if (show_trigger_out && triggers_out_[idx] != 0) {
            glColor3f(0.529f, 0.808f, 0.922f);
            glBegin(GL_LINES);
            glVertex2f(x, -1.0f);
            glVertex2f(x, 1.0f);
            glEnd();
        }
    }

    // Render min/max values using QPainter
    QPainter painter(this);
    painter.setPen(Qt::red);
    QFont font = painter.font();
    font.setPointSize(10);
    painter.setFont(font);

    QString maxText = QString::number(maxVal, 'f', 2);
    QString minText = QString::number(minVal, 'f', 2);

    painter.drawText(width() - 70, 20, maxText);
    painter.drawText(width() - 70, height() - 10, minText);
    
    painter.end();
}

void MainGlWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    update();  // Request a re-draw
}

void MainGlWidget::setSignalSource(int sourceIndex)
{
    if (!dataHandler_) return;

    disconnectCurrentSource();
    
    // Get actual channel counts from dataHandler
    int numRawChannels = dataHandler_->get_channel_count();
    int numROIChannels = dataHandler_->get_ROI_channel_count();

    if (sourceIndex < numRawChannels) currentSignalType = SignalType::EEG;
    else currentSignalType = SignalType::MRI;

    if (sourceIndex >= 0 && sourceIndex < numRawChannels + numROIChannels) {
        // Handle raw channels (0 to numRawChannels-1)
        currentConnection = connect(dataHandler_, &dataHandler::channelDataUpdated,
            this, [this, sourceIndex](int channel, const Eigen::VectorXd &data,
                                    const Eigen::VectorXi &triggersA,
                                    const Eigen::VectorXi &triggersB,
                                    const Eigen::VectorXi &triggersOut,
                                    const Eigen::VectorXd &timeStamps,
                                    const size_t &data_index,
                                    std::string sourceName) {
                if (channel == sourceIndex) {
                    updateData(data, triggersA, triggersB, triggersOut, timeStamps, data_index, sourceName);
                }
            });
    }
}

void MainGlWidget::disconnectCurrentSource()
{
    if (currentConnection) {
        QObject::disconnect(currentConnection);
        currentConnection = QMetaObject::Connection();
    }
}

void MainGlWidget::updateData(const Eigen::VectorXd &newData, 
                            const Eigen::VectorXi &triggers_A, 
                            const Eigen::VectorXi &triggers_B, 
                            const Eigen::VectorXi &triggers_out, 
                            const Eigen::VectorXd &time_stamps,
                            const size_t &data_index,
                            std::string source_name)
{
    if (!pause_view) {
        std::lock_guard<std::mutex> lock(this->dataMutex);
        
        // Check if dimensions differ
        if (dataVector_.size() != newData.size()) dataVector_.resize(newData.size());
        if (triggers_A_.size() != triggers_A.size()) triggers_A_.resize(triggers_A.size());
        if (triggers_B_.size() != triggers_B.size()) triggers_B_.resize(triggers_B.size());
        if (time_stamps_.size() != time_stamps.size()) time_stamps_.resize(time_stamps.size());
        if (triggers_out_.size() != triggers_out.size()) triggers_out_.resize(triggers_out.size());
        
        dataVector_ = newData;
        vectorCapasity = newData.size();
        triggers_A_ = triggers_A;
        triggers_B_ = triggers_B;
        triggers_out_ = triggers_out;
        time_stamps_ = time_stamps;
        data_index_ = data_index;
        source_name_ = QString::fromStdString(source_name);
    }
    update();
}

std::pair<double, double> MainGlWidget::calculateDataRange(const Eigen::VectorXd& data, ScalingMethod method) {
    if (method == ScalingMethod::MINMAX_COEFF) {
        return {data.minCoeff(), data.maxCoeff()};
    }

    // Calculate mean and standard deviation (ignoring zeros)
    double sum = 0;
    int count = 0;
    for (const double& val : data) {
        if (val != 0) {
            sum += val;
            count++;
        }
    }
    double mean = count > 0 ? sum / count : 0;

    // Calculate standard deviation
    double sumSquaredDiff = 0;
    for (const double& val : data) {
        if (val != 0) {
            double diff = val - mean;
            sumSquaredDiff += diff * diff;
        }
    }
    double stdDev = count > 0 ? std::sqrt(sumSquaredDiff / count) : 1.0;

    // Define the range for acceptable values (within 2 standard deviations)
    const double stdDevMultiplier = 2.0;
    double minThreshold = mean - stdDevMultiplier * stdDev;
    double maxThreshold = mean + stdDevMultiplier * stdDev;

    // Find min/max within the acceptable range
    double minVal = std::numeric_limits<double>::max();
    double maxVal = std::numeric_limits<double>::lowest();

    for (const double& val : data) {
        if (val != 0 && val >= minThreshold && val <= maxThreshold) {
            minVal = std::min(minVal, val);
            maxVal = std::max(maxVal, val);
        }
    }

    // Handle edge cases
    if (minVal == std::numeric_limits<double>::max()) {
        minVal = mean - stdDev;
        maxVal = mean + stdDev;
    }

    return {minVal, maxVal};
}