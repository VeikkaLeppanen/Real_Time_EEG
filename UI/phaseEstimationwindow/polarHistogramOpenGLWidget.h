#ifndef POLARHISTOGRAMOPENGLWIDGET_H
#define POLARHISTOGRAMOPENGLWIDGET_H

#include <QOpenGLWidget>
#include <vector>
#include <cmath>
#include <algorithm>
#include <QString>
#include <QFont>
#include <QPainter>

class PolarHistogramOpenGLWidget : public QOpenGLWidget {
    Q_OBJECT

public:
    explicit PolarHistogramOpenGLWidget(QWidget *parent = nullptr)
        : QOpenGLWidget(parent), numBins(36), title1("Pre Initialization"), title2("Post Initialization") {
        bins1.resize(numBins, 0);
        bins2.resize(numBins, 0);
    }

    // Function to add a sample to the first histogram
    void addSampleToFirstCircle(double angle) {
        angle = normalizeAngle(angle);
        int binIndex = static_cast<int>((angle / (2 * M_PI)) * numBins);
        bins1[binIndex] += 1;
        update(); // Request a repaint
    }

    // Function to add a sample to the second histogram
    void addSampleToSecondCircle(double angle) {
        angle = normalizeAngle(angle);
        int binIndex = static_cast<int>((angle / (2 * M_PI)) * numBins);
        bins2[binIndex] += 1;
        update(); // Request a repaint
    }

    // Function to clear samples for both histograms
    void clearSamples() {
        std::fill(bins1.begin(), bins1.end(), 0);
        std::fill(bins2.begin(), bins2.end(), 0);
        update();
    }

    // Function to set titles for both circles
    void setTitles(const QString &t1, const QString &t2) {
        title1 = t1;
        title2 = t2;
        update();
    }

protected:
    void initializeGL() override {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    }

    void resizeGL(int w, int h) override {
        glViewport(0, 0, w, h);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double aspectRatio = static_cast<double>(w) / h;
        if (aspectRatio >= 1.0f) {
            glOrtho(-aspectRatio, aspectRatio, -1.0, 1.0, -1.0, 1.0);
        } else {
            glOrtho(-1.0, 1.0, -1.0 / aspectRatio, 1.0 / aspectRatio, -1.0, 1.0);
        }
        glMatrixMode(GL_MODELVIEW);
    }

    void paintGL() override {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();

        // Draw the first circle and its bars (top)
        drawCircleAndBars(bins1, 0.0, 1.0, 0.9f, 0.6f, 0.8f); // Purple color

        // Draw the second circle and its bars (bottom)
        drawCircleAndBars(bins2, 0.0, -1.0, 0.6f, 0.7f, 0.9f); // Blue color

        // Draw titles using QPainter (must be called outside OpenGL context)
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing);
        painter.setPen(Qt::white);
        painter.setFont(QFont("Helvetica", 12));

        // Draw the title for the first circle (top)
        painter.drawText(width() / 2 - 50, 17, title1);

        // Draw the title for the second circle (bottom)
        painter.drawText(width() / 2 - 50, height() - 5, title2);

        painter.end();
    }

private:
    int numBins; // Number of bins in each histogram
    std::vector<int> bins1; // Counts for the first circle
    std::vector<int> bins2; // Counts for the second circle
    QString title1; // Title for the first circle
    QString title2; // Title for the second circle

    double normalizeAngle(double angle) {
        while (angle < 0) angle += 2 * M_PI;
        return fmod(angle, 2 * M_PI);
    }

    void drawCircleAndBars(const std::vector<int> &bins, double offsetX, double offsetY, float r, float g, float b) {
        double maxCount = *std::max_element(bins.begin(), bins.end());
        if (maxCount == 0) return;

        glPushMatrix();
        glTranslated(offsetX, offsetY, 0);

        // Draw the unit circle
        glColor3f(1.0f, 1.0f, 1.0f);
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < 360; ++i) {
            double angle = i * M_PI / 180.0;
            glVertex2f(cos(angle), sin(angle));
        }
        glEnd();

        // Draw histogram bars
        for (int i = 0; i < numBins; ++i) {
            double angleStart = (i * 2 * M_PI) / numBins;
            double angleEnd = ((i + 1) * 2 * M_PI) / numBins;
            double normalizedHeight = static_cast<double>(bins[i]) / maxCount;

            glColor3f(r, g, b);
            glBegin(GL_TRIANGLE_FAN);
            glVertex2f(0.0f, 0.0f);
            for (double angle = angleStart; angle <= angleEnd; angle += 0.01) {
                double x = normalizedHeight * cos(angle);
                double y = normalizedHeight * sin(angle);
                glVertex2f(x, y);
            }
            glVertex2f(normalizedHeight * cos(angleEnd), normalizedHeight * sin(angleEnd));
            glEnd();

            // Draw the outline of the bar
            glColor3f(1.0f, 1.0f, 1.0f);
            glBegin(GL_LINE_STRIP);
            glVertex2f(0.0f, 0.0f);
            for (double angle = angleStart; angle <= angleEnd; angle += 0.01) {
                double x = normalizedHeight * cos(angle);
                double y = normalizedHeight * sin(angle);
                glVertex2f(x, y);
            }
            glVertex2f(normalizedHeight * cos(angleEnd), normalizedHeight * sin(angleEnd));
            glVertex2f(0.0f, 0.0f);
            glEnd();
        }

        glPopMatrix();
    }
};

#endif // POLARHISTOGRAMOPENGLWIDGET_H
