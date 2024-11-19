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

enum EditorMode {
  BRUSH,
  ERASE,
  RECTANGLE
};

enum SliceType {
    AXIAL,
    CORONAL,
    SAGITTAL
};

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

signals:
    void ROI_update(const std::vector<std::string> &names, std::vector<bool> ROI_toggleStatus, std::vector<bool> ROI_visibility);

public slots:
    void addButton_clicked();
    void loadButton_clicked(const QString& filePath);
    void saveButton_clicked();
    void deleteButton_clicked();
    void visibleButton_clicked();

    void updateToggleStates(std::vector<bool> states) { ROI_toggle_states = states; };

    void onSliderValueChanged(int value);
    void colorButton_clicked(QColor color);
    
    void editButton_toggled(bool checked);
    void undoButton_clicked();
    void redoButton_clicked();

    void brushButton_clicked() { editorMode = BRUSH; }
    void eraseButton_clicked() { editorMode = ERASE; }
    void rectangleButton_clicked() { editorMode = RECTANGLE; }

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
    std::vector<bool> ROI_toggle_states;
    std::vector<bool> ROI_visibility;
    std::vector<QColor> ROI_colors;
    std::vector<float> ROI_opacities;

    // Slice indices for each plane T1
    int i; // Sagittal (along x-axis)
    int j; // Coronal (along y-axis)
    int k; // Axial (along z-axis)

    std::vector<std::string> T1_orientation = {"R", "A", "S"};
    // std::vector<std::string> T1_orientation = {"L", "P", "I"};
    std::vector<float> T1_pixDims = {1.0, 1.0, 1.0};

    // Member variables to store projection parameters
    float leftAxial, rightAxial, bottomAxial, topAxial;
    float leftCoronal, rightCoronal, bottomCoronal, topCoronal;
    float leftSagittal, rightSagittal, bottomSagittal, topSagittal;

    float imgWidth_axial;
    float imgHeight_axial;
    float imgWidth_coronal;
    float imgHeight_coronal;
    float imgWidth_sagittal;
    float imgHeight_sagittal;

    Eigen::Matrix4f constructMatrix(float ijk2xyz[3][4]);

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

    // ROI parameters
    bool editMode = false;
    EditorMode editorMode = BRUSH;
    QPoint rectStartPoint;   // Starting point of the rectangle
    QPoint rectCurrentPoint; // Current point as the mouse moves
    QPoint rectStartImageCoords;
    QPoint rectEndImageCoords;
    SliceType currentSliceType;
    bool isDrawingRect = false; // Flag to indicate if rectangle drawing is in progress
    void applyRectangleToROI(const QPoint& startImageCoords, const QPoint& endImageCoords);
    QPoint screenToImageCoordinates(const QPoint& screenPoint);
    void convertAxialScreenToImage(const QPoint& screenPoint, int viewportWidth, int viewportHeight, int& imgi, int& imgj);
    void convertCoronalScreenToImage(const QPoint& screenPoint, int viewportWidth, int viewportHeight, int& imgi, int& imgk);
    void convertSagittalScreenToImage(const QPoint& screenPoint, int viewportWidth, int viewportHeight, int& imgj, int& imgk);
};

#endif // MAINGLWIDGET_H
