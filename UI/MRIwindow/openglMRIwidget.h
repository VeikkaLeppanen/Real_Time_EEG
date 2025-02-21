#ifndef OPENGLMRIWIDGET_H
#define OPENGLMRIWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QTimer>
#include <iostream>
#include <float.h>
#include <unordered_set>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QPainter>
#include <QFont>
#include <QMessageBox>
#include <image/image.h>
#include <image/image_operators.h>
#include <image/orientation.h>
#include <Eigen/Dense>
#include <QFileSystemWatcher>
#include <QDir>
#include <QStringList>
#include "../../dataHandler/dataHandler.h"
#include <QtConcurrent>
#include <QFutureWatcher>

enum BrushMode {
  BRUSH_SQUARE,
  BRUSH_CIRCLE,
  RECTANGLE
};

enum SliceType {
    AXIAL,
    CORONAL,
    SAGITTAL
};

struct Viewport {
    int x;
    int y;
    int width;
    int height;
};

class OpenGLMRIWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit OpenGLMRIWidget(QWidget *parent = nullptr);
    void setSliceIndices(int new_i, int new_j, int new_k);
    void loadImage_T1(const QString& filePath);
    void loadImage_fMRI(const QString& filePath);
    void setWatchFolder(const QString& path);
    void setDataHandler(dataHandler* handler) { data_handler = handler; }
    const std::vector<std::string>& getROINames() const { return ROI_names; }

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void drawDirectionLabels(QPainter &painter, const Viewport &viewport, const QString &sliceType);
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

signals:
    void ROI_update(const std::vector<std::string> &names, std::vector<bool> ROI_toggleStatus);
    void ROIListChanged(const std::vector<std::string>& channelNames);

public slots:
    void addButton_clicked();
    void loadButton_clicked(const QString& filePath);
    void saveButton_clicked();
    void deleteButton_clicked();

    void updateToggleStates(std::vector<bool> states) { ROI_toggle_states = states; update(); };
    void updateTargetROI(int index) { target_ROI = index; };

    void onSliderValueChanged(int value);
    void colorButton_clicked(QColor color);
    
    void editButton_toggled(bool checked);
    void undoButton_clicked();
    void redoButton_clicked();

    void onMarkButtonClicked() { isMarking = true; }
    void onEraseButtonClicked() { isMarking = false; }
    void onBrushSelected(int index) {
        if (index >= 0 && index <= RECTANGLE) {
            editorMode = static_cast<BrushMode>(index);
        } else {
            qWarning() << "Invalid brush mode index selected:" << index;
        }
    }
    void onBrushSizeChanged(int size) { brush_size = size; }

    void aboveButton_clicked();
    void belowButton_clicked();

    void edit_ROIs();

    void reset_undo_stacks() {
        undo_stack_temp.clear();
        undo_stack.clear();
        redo_stack.clear();
    }

private slots:
    void updateGraph();
    void handleNewFile(const QString& path);

private:
    const bool debug = false;
    void print_debug(std::string msg) {
        if (debug) std::cout << msg << std::endl;
    };
    
    // MRI image data
    NIBR::Image<float> T1_image;
    NIBR::Image<float> fMRI_image;
    NIBR::Image<float> fMRI_image_resampled;
    float minValue;
    float maxValue;
    float mean_intensity;
    float minValue_f;
    float maxValue_f;
    float mean_intensity_f;

    float T1_contrast = 1.0f;
    float fMRI_contrast = 1.0f;
    float fMRI_overlayOpacity; // Value between 0.0f and 1.0f
    QPoint lastContrastPoint;

    int target_ROI;
    std::vector<std::string> ROI_names;
    std::vector<std::shared_ptr<NIBR::Image<bool>>> ROI_vector;
    std::vector<bool> ROI_toggle_states;
    std::vector<QColor> ROI_colors;
    std::vector<float> ROI_opacities;

    std::unordered_set<int64_t> undo_stack_temp;
    std::vector<std::vector<int64_t>> undo_stack;
    std::vector<std::vector<int64_t>> redo_stack;
    void revert_ROI_at_index(int64_t image_index);
    void paint_square_ROI(int brush_size);
    void paint_circle_ROI(int brush_size);
    void set_ROI_at_index(int64_t image_index, bool ROI_value);

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
    BrushMode editorMode = BRUSH_SQUARE;
    int brush_size = 0;
    bool isMarking = true;
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

    // Add these new member variables
    Eigen::VectorXd ROI_fMRI_means;  // Buffer to store mean fMRI values for each ROI
    void calculateROIMeans();  // Helper function to calculate means

    QFileSystemWatcher* fileWatcher;
    QString watchFolderPath;
    QStringList processedFiles;  // Keep track of processed files

    dataHandler* data_handler = nullptr;
};

#endif // OPENGLMRIWIDGET_H
