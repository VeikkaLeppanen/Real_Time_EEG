#include "glwidget.h"

Glwidget::Glwidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
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
    glClear(GL_COLOR_BUFFER_BIT);
    // Add your OpenGL drawing code here to visualize the EEG data
    // For example, drawing a simple line graph or points

    // Example to plot random data (replace with actual data handling)
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0, 1.0, 1.0); // White color for the lines
    for (int i = 0; i < 100; ++i) {
        glVertex2f(i / 50.0f - 1.0f, (rand() % 100) / 50.0f - 1.0f);
    }
    glEnd();
}

void Glwidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    update();  // Request a re-draw
}
