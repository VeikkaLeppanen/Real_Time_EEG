#include "glwidget.h"

Glwidget::Glwidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    QTimer *timer = new QTimer(this);
    newMatrix_ = Eigen::MatrixXd::Zero(1, 10000);
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
    if (channelCheckStates_.size() == 0) return;  // Ensure there is data to draw

    int windowHeight = height();
    int rowHeight = windowHeight / std::count(channelCheckStates_.begin(), channelCheckStates_.end(), true);  // Divide the window height by the number of rows
    // std::cout << std::count(channelCheckStates_.begin(), channelCheckStates_.end(), true) << '\n';
    glClear(GL_COLOR_BUFFER_BIT);

    int graph_index = 0;
    for (int row = 0; row < newMatrix_.rows(); row++) {
        if (!channelCheckStates_[row]) continue;
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
        Eigen::VectorXd dataVector = newMatrix_.row(n_channels - 1 - row);
        double minVal = dataVector.minCoeff();
        double maxVal = dataVector.maxCoeff();

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

    // Overlay text using QPainter
    QPainter painter(this);
    painter.setPen(Qt::red);
    painter.setFont(QFont("Arial", 10)); // Set font here

    graph_index = 0;
    for (int row = 0; row < newMatrix_.rows(); row++) {
        if (!channelCheckStates_[row]) continue;
        // Calculate position for text
        int yPos = windowHeight - (rowHeight * graph_index + rowHeight / 2 + 10);  // Adjust vertical position
        QString name = "Undefined";  // Default name if no channel name is available
        if (row < channelNames_.size()) {
            name = channelNames_.at(n_channels - 1 - row);
        }
        painter.drawText(10, yPos, name);
        graph_index++;
    }

    painter.end(); // Ensure QPainter is properly closed
}

void Glwidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    emit fetchData();
    update();  // Request a re-draw
}
