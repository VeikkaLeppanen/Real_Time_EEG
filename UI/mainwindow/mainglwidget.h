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
#include <image/image_operators.h>
#include <image/orientation.h>
#include <Eigen/Dense>

class MainGlWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit MainGlWidget(QWidget *parent = nullptr);
    void setSliceIndices(int new_i, int new_j, int new_k);
    void loadImage_T1(const QString& filePath);
    void loadImage_fMRI(const QString& filePath);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

    // void resizeEvent(QResizeEvent *event) override {
    //     int newWidth = event->size().width();
    //     int newHeight = event->size().height();

    //     // Calculate the maximum size that fits within the new dimensions
    //     if (newWidth > newHeight * aspectRatio) {
    //         // Height is the limiting dimension
    //         newWidth = static_cast<int>(newHeight * aspectRatio);
    //     } else {
    //         // Width is the limiting dimension
    //         newHeight = static_cast<int>(newWidth / aspectRatio);
    //     }

    //     // Set the viewport and resize the widget
    //     setGeometry((width() - newWidth) / 2, (height() - newHeight) / 2, newWidth, newHeight);
    //     QOpenGLWidget::resizeEvent(event);
    // }

signals:
    void ROI_update(const std::vector<std::string> &names, const std::vector<bool> &ROI_visibility);

public slots:
    void addButton_clicked();
    void loadButton_clicked(const QString& filePath);
    void saveButton_clicked(const std::vector<bool>& states);
    void deleteButton_clicked(const std::vector<bool>& states);
    void visibleButton_clicked(const std::vector<bool>& states);

private slots:
    void updateGraph();

private:
    // MRI image data
    NIBR::Image<float> T1_image;
    NIBR::Image<float> fMRI_image;
    NIBR::Image<float> fMRI_image_resampled;
    float minValue;
    float maxValue;
    float minValue_f;
    float maxValue_f;

    std::vector<std::string> ROI_names;
    std::vector<NIBR::Image<bool>> ROI_vector;
    std::vector<bool> ROI_visibility;

    // Slice indices for each plane T1
    int i; // Sagittal (along x-axis)
    int j; // Coronal (along y-axis)
    int k; // Axial (along z-axis)

    std::vector<std::string> T1_orientation = {"R", "A", "S"};
    // std::vector<std::string> T1_orientation = {"L", "P", "I"};
    std::vector<float> T1_pixDims = {1.0, 1.0, 1.0};
    // std::vector<float> T1_pixDims = {2.0, 1.5, 0.5};

    Eigen::Matrix4f constructMatrix(float ijk2xyz[3][4]);
    float getInterpolatedVoxelValue(float* data, float x, float y, float z, int dimX, int dimY, int dimZ);

    // Zoom and Pan
    float zoomFactor;
    QVector3D panOffset;
    QPoint lastPanPoint; // For tracking mouse movement during panning
    float overlayOpacity; // Value between 0.0f and 1.0f
    // const double aspectRatio; // Set your aspect ratio (e.g., 16:9 or 4:3)

    // Mouse event handlers
    void handleAxialClick(const QPoint& mousePos, int viewportWidth, int viewportHeight);
    void handleCoronalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight);
    void handleSagittalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight);

    // Scroll event handlers
    void handleAxialScroll(QWheelEvent *event);
    void handleCoronalScroll(QWheelEvent *event);
    void handleSagittalScroll(QWheelEvent *event);

    // Mouse state
    bool mousePressed;
    Qt::MouseButton pressedButton;
};

#endif // MAINGLWIDGET_H
