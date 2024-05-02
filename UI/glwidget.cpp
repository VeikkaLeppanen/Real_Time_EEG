#include "glwidget.h"

Glwidget::Glwidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    QTimer *timer = new QTimer(this);
    dataMatrix = Eigen::MatrixXd::Zero(10, 10000);
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
    if (dataMatrix.rows() == 0) return;  // Ensure there is data to draw

    int windowHeight = height();
    int rowHeight = windowHeight / dataMatrix.rows();  // Divide the window height by the number of rows

    glClear(GL_COLOR_BUFFER_BIT);

    for (int row = 0; row < dataMatrix.rows(); ++row) {
        // Set viewport for this row
        glViewport(0, row * rowHeight, width(), rowHeight);

        // Set up the projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);  // Set the coordinate system to cover [-1,1] in both axes

        // Switch back to modelview matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Prepare data and draw the line strip
        Eigen::VectorXd dataVector = dataMatrix.row(row);
        double minVal = dataVector.minCoeff();
        double maxVal = dataVector.maxCoeff();

        glColor3f(1.0, 1.0, 1.0); // Set the color to white for the lines
        glBegin(GL_LINE_STRIP);

        for (int i = 0; i < dataVector.size(); ++i) {
            float x = (float)i / (dataVector.size() - 1) * 2.0f - 1.0f;  // Scale x to [-1,1]
            float y = 0;
            if (maxVal != minVal) {  // Avoid division by zero
                y = (dataVector(i) - minVal) / (maxVal - minVal) * 2.0f - 1.0f;  // Normalize and scale y to [-1,1]
            }
            glVertex2f(x, y);
        }
        glEnd();
    }
}



void Glwidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    emit fetchData();
    update();  // Request a re-draw
}
