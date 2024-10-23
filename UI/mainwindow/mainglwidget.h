#ifndef MAINGLWIDGET_H
#define MAINGLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QTimer>
#include <iostream>
#include <float.h>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QPainter>
#include <QFont>
#include <image/image.h>

class MainGlWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit MainGlWidget(QWidget *parent = nullptr);
    void setSliceIndices(int new_i, int new_j, int new_k);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override; // Add this line
    void mouseReleaseEvent(QMouseEvent *event) override; // Add this line
    void wheelEvent(QWheelEvent *event) override; // Add this line
    void keyPressEvent(QKeyEvent *event) override;

private slots:
    void updateGraph();

private:
    // MRI image data
    NIBR::Image<float> dMRI_image;
    float minValue;
    float maxValue;

    // Slice indices for each plane
    int i; // Sagittal (along x-axis)
    int j; // Coronal (along y-axis)
    int k; // Axial (along z-axis)
    
    // Zoom and Pan
    float zoomFactor;
    QVector3D panOffset;
    QPoint lastPanPoint; // For tracking mouse movement during panning

    // Mouse event handlers
    void handleAxialClick(const QPoint& mousePos, int viewportWidth, int viewportHeight);
    void handleCoronalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight);
    void handleSagittalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight);

    // Scroll event handlers
    void handleAxialScroll(QWheelEvent *event);
    void handleCoronalScroll(QWheelEvent *event);
    void handleSagittalScroll(QWheelEvent *event);
    
    // Mouse state
    bool mousePressed; // Add this line
    Qt::MouseButton pressedButton; // Add this line
};

#endif // MAINGLWIDGET_H
