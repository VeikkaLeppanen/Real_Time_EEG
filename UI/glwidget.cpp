#include "glwidget.h"

Glwidget::Glwidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    windowLength_seconds = 2;   // Total time in seconds displayed
    time_line_spacing = 500;    // Line spacing in milliseconds
    totalTimeLines = (windowLength_seconds * 1000) / time_line_spacing;  // Calculate number of lines to draw
    QTimer *timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &Glwidget::updateGraph);
    timer->start(16); // Update approximately every 16 ms (60 FPS)
}

void Glwidget::initializeGL()
{
    initializeOpenGLFunctions();
    glClearColor(0.0, 0.0, 0.0, 1.0); // Set the clearing color to black
}

void Glwidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h); // Set the viewport to cover the entire widget
}

void Glwidget::paintGL()
{
    if (dataMatrix_.rows() == 0) return;

    // Initializing positional parameters
    int windowHeight = height();
    int rowHeight = windowHeight / n_channels;

    // min max values for y axis
    Eigen::VectorXd min_coeffs = Eigen::VectorXd::Zero(n_channels);
    Eigen::VectorXd max_coeffs = Eigen::VectorXd::Zero(n_channels);
    
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
        Eigen::VectorXd dataVector = dataMatrix_.row(n_channels - 1 - row);
        double minVal = dataVector.minCoeff();
        min_coeffs(row) = minVal;
        double maxVal = dataVector.maxCoeff();
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
                float x = (float)i / (totalDataPoints - 1) * 2.0f - 1.0f;
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
                float x = (float)i / (totalDataPoints - 1) * 2.0f - 1.0f;
                glColor3f(0.0, 1.0, 0.0);
                glBegin(GL_LINES);
                glVertex2f(x, -1.0);
                glVertex2f(x, 1.0);
                glEnd();
            }
        }
    }

    if (drawXaxis) {
        // Enable stipple for dashed lines
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1, 0x00FF);  // 1x repeat factor, 0x00FF pattern
        glColor3f(0.5, 0.5, 0.5);  // Gray color for the lines

        // Draw each vertical line
        for (int i = 0; i < totalTimeLines; i++) {
            float x = ((float)i / totalTimeLines) * 2.0f - 1.0f;  // Convert index to OpenGL coordinates

            // Check if the current line is at a whole second
            if ((i * time_line_spacing) % 1000 == 0) {
                glDisable(GL_LINE_STIPPLE);  // Disable stipple for whole second lines
                glColor3f(1.0, 0.0, 0.0);    // Red color for whole second lines
            } else {
                glEnable(GL_LINE_STIPPLE);
                glLineStipple(1, 0x00FF);   // 1x repeat factor, 0x00FF pattern
                glColor3f(0.5, 0.5, 0.5);   // Gray color for non-whole second lines
            }

            glBegin(GL_LINES);
            glVertex2f(x, -1.0);
            glVertex2f(x, 1.0);
            glEnd();
        }

        glDisable(GL_LINE_STIPPLE); // Disable stipple after drawing lines
    }

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
            name = channelNames_.at(n_channels - 1 - row);
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

void Glwidget::updateMatrix(const Eigen::MatrixXd &newMatrix, const Eigen::VectorXi &triggers_A, const Eigen::VectorXi &triggers_B) { 
    if (!pause_view) {
        dataMatrix_ = newMatrix;
        matrixCapasity = newMatrix.cols();
        n_channels = newMatrix.rows();
        triggers_A_ = triggers_A;
        triggers_B_ = triggers_B;
    }
}

void Glwidget::updateChannelDisplayState(std::vector<bool> channelCheckStates) {
    channelCheckStates_ = channelCheckStates; 
}

void Glwidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    emit fetchData();
    update();  // Request a re-draw
}

void Glwidget::updateChannelNamesQt(QStringList channelNames) {
    channelNames_ = channelNames; 
}

void Glwidget::updateChannelNamesSTD(std::vector<std::string> channelNames) {
    QStringList newNames;
    for(size_t i = 0; i < channelNames.size(); i++) {
        newNames.append(QString::fromStdString(channelNames[i])); 
    }
    channelNames_ = newNames;
}
