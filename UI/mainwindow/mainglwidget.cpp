#include "mainglwidget.h"

MainGlWidget::MainGlWidget(QWidget *parent)
    : QOpenGLWidget(parent), 
      dMRI_image("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/devices/MRI/T1.nii.gz"), 
      mousePressed(false),
      zoomFactor(1.0f),
      panOffset(0.0f, 0.0f, 0.0f) // Initialize panOffset as QVector3D
{
    QTimer *timer = new QTimer(this);

    std::cout << "Loading dMRI image" << std::endl;
    dMRI_image.read();
    dMRI_image.printInfo();

    // Compute min and max voxel values
    minValue = FLT_MAX;
    maxValue = -FLT_MAX;
    for (int idx = 0; idx < dMRI_image.voxCnt; ++idx) {
        float value = dMRI_image.data[idx];
        if (value < minValue) minValue = value;
        if (value > maxValue) maxValue = value;
    }
    std::cout << "Voxel value range: [" << minValue << ", " << maxValue << "]" << std::endl;

    // Set default slice indices to the middle of the image
    i = dMRI_image.imgDims[0] / 2; // Sagittal slice index
    j = dMRI_image.imgDims[1] / 2; // Coronal slice index
    k = dMRI_image.imgDims[2] / 2; // Axial slice index
    
    this->setFocusPolicy ( Qt::StrongFocus );
    
    connect(timer, &QTimer::timeout, this, &MainGlWidget::updateGraph);
    timer->start(16); // Update approximately every 16 ms (60 FPS)
}

void MainGlWidget::setSliceIndices(int new_i, int new_j, int new_k)
{
    i = new_i;
    j = new_j;
    k = new_k;
    update(); // Trigger a repaint
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Dimensions of the MRI image
    int dimX = dMRI_image.imgDims[0];
    int dimY = dMRI_image.imgDims[1];
    int dimZ = dMRI_image.imgDims[2];

    float* voxelData = dMRI_image.data;

    // Function to compute the 1D index from 3D coordinates
    auto voxelIndex = [&](int x, int y, int z) -> int {
        return x + y * dimX + z * dimX * dimY;
    };

    // Width and height of the widget
    int widgetWidth = width();
    int widgetHeight = height();

    // Calculate the width of each viewport section
    int viewportWidth = widgetWidth / 3;

    // Define viewports for each slice in the order: Axial, Coronal, Sagittal
    struct Viewport {
        int x;
        int y;
        int width;
        int height;
    };

    Viewport axialViewport = {0, 0, viewportWidth, widgetHeight};
    Viewport coronalViewport = {viewportWidth, 0, viewportWidth, widgetHeight};
    Viewport sagittalViewport = {2 * viewportWidth, 0, viewportWidth, widgetHeight};

    // Helper function to set up projection and modelview matrices
    auto setupOrtho = [&](int imgWidth, int imgHeight, float panX, float panY) {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
    
        float halfWidth = imgWidth / 2.0f;
        float halfHeight = imgHeight / 2.0f;
    
        // Calculate the zoomed and panned view
        float left = (-halfWidth - panX) / zoomFactor;
        float right = (halfWidth - panX) / zoomFactor;
        float bottom = (-halfHeight - panY) / zoomFactor;
        float top = (halfHeight - panY) / zoomFactor;
    
        glOrtho(left, right, bottom, top, -1, 1);
    
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    };

    // --------------------
    // Draw Axial Slice
    // --------------------
    glViewport(axialViewport.x, axialViewport.y, axialViewport.width, axialViewport.height);

    // Use the panOffset and zoomFactor
    setupOrtho(dimX, dimY, panOffset.x(), panOffset.y());

    int z = std::clamp(k, 0, dimZ - 1);

    // Draw the axial slice
    glBegin(GL_QUADS);
    for (int y = 0; y < dimY; ++y) {
        int flippedY = dimY - y - 1; // Flip the y-coordinate
        for (int x = 0; x < dimX; ++x) {
            int idx = voxelIndex(x, flippedY, z);
            float value = voxelData[idx];

            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            glColor3f(gray, gray, gray);

            // Adjust coordinates to center around (0,0)
            float x0 = (x - dimX / 2.0f) + 0.5f;
            float y0 = (y - dimY / 2.0f) + 0.5f;
            float x1 = x0 + 1.0f;
            float y1 = y0 + 1.0f;

            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();

    // Draw lines on axial slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    // Vertical line at x = i
    float lineX = (i - dimX / 2.0f) + 1.0f;
    glVertex2f(lineX, -dimY / 2.0f);
    glVertex2f(lineX, dimY / 2.0f);
    // Horizontal line at y = j
    float lineY = ((dimY - j - 1) - dimY / 2.0f) + 1.0f;
    glVertex2f(-dimX / 2.0f, lineY);
    glVertex2f(dimX / 2.0f, lineY);
    glEnd();

    // --------------------
    // Draw Coronal Slice
    // --------------------
    glViewport(coronalViewport.x, coronalViewport.y, coronalViewport.width, coronalViewport.height);

    // Use the panOffset and zoomFactor
    setupOrtho(dimX, dimZ, panOffset.x(), panOffset.z());

    int y = std::clamp(j, 0, dimY - 1);

    // Draw the coronal slice
    glBegin(GL_QUADS);
    for (int z = 0; z < dimZ; ++z) {
        int flippedZ = dimZ - z - 1; // Flip the z-coordinate
        for (int x = 0; x < dimX; ++x) {
            int idx = voxelIndex(x, y, z);
            float value = voxelData[idx];

            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            glColor3f(gray, gray, gray);

            // Adjust coordinates to center around (0,0)
            float x0 = (x - dimX / 2.0f) + 0.5f;
            float y0 = (z - dimZ / 2.0f) + 0.5f;
            float x1 = x0 + 1.0f;
            float y1 = y0 + 1.0f;

            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();

    // Draw lines on coronal slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    // Vertical line at x = i
    lineX = (i - dimX / 2.0f) + 1.0f;
    glVertex2f(lineX, -dimZ / 2.0f);
    glVertex2f(lineX, dimZ / 2.0f);
    // Horizontal line at z = k
    lineY = ((k - dimZ / 2.0f) + 1.0f);
    glVertex2f(-dimX / 2.0f, lineY);
    glVertex2f(dimX / 2.0f, lineY);
    glEnd();

    // --------------------
    // Draw Sagittal Slice
    // --------------------
    glViewport(sagittalViewport.x, sagittalViewport.y, sagittalViewport.width, sagittalViewport.height);

    // Use the panOffset and zoomFactor
    setupOrtho(dimY, dimZ, panOffset.y(), panOffset.z());

    int x = std::clamp(i, 0, dimX - 1);

    // Draw the sagittal slice
    glBegin(GL_QUADS);
    for (int z = 0; z < dimZ; ++z) {
        for (int y = 0; y < dimY; ++y) {
            int idx = voxelIndex(x, y, z);
            float value = voxelData[idx];

            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            glColor3f(gray, gray, gray);

            // Adjust coordinates to center around (0,0)
            float x0 = (y - dimY / 2.0f) + 0.5f;
            float y0 = (z - dimZ / 2.0f) + 0.5f;
            float x1 = x0 + 1.0f;
            float y1 = y0 + 1.0f;

            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();

    // Draw lines on sagittal slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    // Vertical line at y = j
    lineX = (j - dimY / 2.0f) + 1.0f;
    glVertex2f(lineX, -dimZ / 2.0f);
    glVertex2f(lineX, dimZ / 2.0f);
    // Horizontal line at z = k
    lineY = ((k - dimZ / 2.0f) + 1.0f);
    glVertex2f(-dimY / 2.0f, lineY);
    glVertex2f(dimY / 2.0f, lineY);
    glEnd();

    // Retrieve the voxel value
    int ix = std::clamp(i, 0, dimX - 1);
    int iy = std::clamp(j, 0, dimY - 1);
    int iz = std::clamp(k, 0, dimZ - 1);

    // Retrieve the target voxel value
    int trgt_idx = voxelIndex(ix, iy, iz);
    float trgt_voxelValue = voxelData[trgt_idx];

    // Begin QPainter
    QPainter painter(this);

    // Set text properties
    painter.setPen(Qt::yellow);
    painter.setFont(QFont("Arial", 12));

    // Prepare the text
    QString text = QString("Value: %1\nVoxel: (%2, %3, %4)")
        .arg(trgt_voxelValue).arg(ix).arg(iy).arg(iz);

    // Position the text at the bottom left corner
    QRect textRect(10, height() - 60, 200, 50); // x, y, width, height

    painter.drawText(textRect, Qt::AlignLeft | Qt::AlignBottom | Qt::TextWordWrap, text);

    // End QPainter
    painter.end();
}

void MainGlWidget::handleAxialClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Dimensions of the Axial slice
    int dimX = dMRI_image.imgDims[0];
    int dimY = dMRI_image.imgDims[1];

    // Calculate the position within the axial viewport
    int xInViewport = mousePos.x();
    int yInViewport = viewportHeight - mousePos.y(); // Flip y-coordinate

    // Convert viewport coordinates to Normalized Device Coordinates (NDC)
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewport) / viewportHeight) * 2.0f - 1.0f;

    // Undo the projection transformation (glOrtho)
    float halfWidth = dimX / 2.0f;
    float halfHeight = dimY / 2.0f;

    float panX = panOffset.x();
    float panY = panOffset.y();

    float left = (-halfWidth - panX) / zoomFactor;
    float right = (halfWidth - panX) / zoomFactor;
    float bottom = (-halfHeight - panY) / zoomFactor;
    float top = (halfHeight - panY) / zoomFactor;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (right - left) + left;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (top - bottom) + bottom;

    // Convert world coordinates to image coordinates
    int imgX = static_cast<int>(worldX + dimX / 2.0f - 0.5f);
    int imgY = static_cast<int>(worldY + dimY / 2.0f - 0.5f);

    // Flip imgY due to image flipping
    imgY = dimY - imgY - 1;

    // Clamp the indices to valid ranges
    i = std::clamp(imgX, 0, dimX - 1);
    j = std::clamp(imgY, 0, dimY - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::handleCoronalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Dimensions of the Coronal slice
    int dimX = dMRI_image.imgDims[0];
    int dimZ = dMRI_image.imgDims[2];

    // Adjust mouse position for the coronal viewport
    int xInViewport = mousePos.x() - viewportWidth; // Subtract viewport offset
    int yInViewport = viewportHeight - mousePos.y(); // Flip y-coordinate

    // Convert viewport coordinates to NDC
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewport) / viewportHeight) * 2.0f - 1.0f;

    // Undo the projection transformation
    float halfWidth = dimX / 2.0f;
    float halfHeight = dimZ / 2.0f;

    float panX = panOffset.x();
    float panZ = panOffset.z();

    float left = (-halfWidth - panX) / zoomFactor;
    float right = (halfWidth - panX) / zoomFactor;
    float bottom = (-halfHeight - panZ) / zoomFactor;
    float top = (halfHeight - panZ) / zoomFactor;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (right - left) + left;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (top - bottom) + bottom;

    // Convert world coordinates to image coordinates
    int imgX = static_cast<int>(worldX + dimX / 2.0f - 0.5f);
    int imgZ = static_cast<int>(worldY + dimZ / 2.0f - 0.5f);

    // Clamp the indices to valid ranges
    i = std::clamp(imgX, 0, dimX - 1);
    k = std::clamp(imgZ, 0, dimZ - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::handleSagittalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Dimensions of the Sagittal slice
    int dimY = dMRI_image.imgDims[1];
    int dimZ = dMRI_image.imgDims[2];

    // Adjust mouse position for the sagittal viewport
    int xInViewport = mousePos.x() - 2 * viewportWidth; // Subtract viewport offset
    int yInViewport = viewportHeight - mousePos.y(); // Flip y-coordinate

    // Convert viewport coordinates to NDC
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewport) / viewportHeight) * 2.0f - 1.0f;

    // Undo the projection transformation
    float halfWidth = dimY / 2.0f;
    float halfHeight = dimZ / 2.0f;

    float panY = panOffset.y();
    float panZ = panOffset.z();

    float left = (-halfWidth - panY) / zoomFactor;
    float right = (halfWidth - panY) / zoomFactor;
    float bottom = (-halfHeight - panZ) / zoomFactor;
    float top = (halfHeight - panZ) / zoomFactor;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (right - left) + left;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (top - bottom) + bottom;

    // Convert world coordinates to image coordinates
    int imgY = static_cast<int>(worldX + dimY / 2.0f - 0.5f);
    int imgZ = static_cast<int>(worldY + dimZ / 2.0f - 0.5f);

    // Clamp the indices to valid ranges
    j = std::clamp(imgY, 0, dimY - 1);
    k = std::clamp(imgZ, 0, dimZ - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::mousePressEvent(QMouseEvent *event)
{
    mousePressed = true;
    pressedButton = event->button();

    if (pressedButton == Qt::RightButton || event->modifiers() & Qt::ControlModifier)
    {
        // Start panning
        lastPanPoint = event->pos();
    }
    else if (pressedButton == Qt::LeftButton)
    {
        // Handle crosshair movement
        QPoint mousePos = event->pos();
        int widgetWidth = width();
        int viewportWidth = widgetWidth / 3;

        if (mousePos.x() < viewportWidth) {
            handleAxialClick(mousePos, viewportWidth, height());
        } else if (mousePos.x() < 2 * viewportWidth) {
            handleCoronalClick(mousePos, viewportWidth, height());
        } else {
            handleSagittalClick(mousePos, viewportWidth, height());
        }
    }
}

void MainGlWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (!mousePressed) {
        return;
    }

    if (pressedButton == Qt::RightButton || event->modifiers() & Qt::ControlModifier) {
        // Handle panning
        QPoint currentPoint = event->pos();
        QPoint delta = -(currentPoint - lastPanPoint);
        lastPanPoint = currentPoint;

        // Determine which viewport the cursor is over
        int widgetWidth = width();
        int widgetHeight = height();
        int viewportWidth = widgetWidth / 3;
        int viewportIndex = currentPoint.x() / viewportWidth;

        float deltaX = delta.x();
        float deltaY = -delta.y(); // Invert y to match OpenGL coordinate system

        float scaleX, scaleY;

        if (viewportIndex == 0) {
            // Axial slice (x, y plane)
            scaleX = (dMRI_image.imgDims[0] + zoomFactor) / viewportWidth;
            scaleY = (dMRI_image.imgDims[1] + zoomFactor) / widgetHeight;

            panOffset.setX(panOffset.x() - deltaX * scaleX);
            panOffset.setY(panOffset.y() - deltaY * scaleY);

        } else if (viewportIndex == 1) {
            // Coronal slice (x, z plane)
            scaleX = (dMRI_image.imgDims[0] + zoomFactor) / viewportWidth;
            scaleY = (dMRI_image.imgDims[2] + zoomFactor) / widgetHeight;

            panOffset.setX(panOffset.x() - deltaX * scaleX);
            panOffset.setZ(panOffset.z() - deltaY * scaleY);

        } else {
            // Sagittal slice (y, z plane)
            scaleX = (dMRI_image.imgDims[1] + zoomFactor) / viewportWidth;
            scaleY = (dMRI_image.imgDims[2] + zoomFactor) / widgetHeight;

            panOffset.setY(panOffset.y() - deltaX * scaleX);
            panOffset.setZ(panOffset.z() - deltaY * scaleY);
        }

        update();

    } else if (pressedButton == Qt::LeftButton) {
        // Handle crosshair movement
        QPoint mousePos = event->pos();
        int widgetWidth = width();
        int viewportWidth = widgetWidth / 3;

        if (mousePos.x() < viewportWidth) {
            handleAxialClick(mousePos, viewportWidth, height());
        } else if (mousePos.x() < 2 * viewportWidth) {
            handleCoronalClick(mousePos, viewportWidth, height());
        } else {
            handleSagittalClick(mousePos, viewportWidth, height());
        }
    }
}

void MainGlWidget::mouseReleaseEvent(QMouseEvent *event)
{
    Q_UNUSED(event);
    // Reset the mouse pressed flag
    mousePressed = false;
}

void MainGlWidget::handleAxialScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120; // Each notch is 15 degrees; 120 units per notch
    k = std::clamp(k + delta, 0, static_cast<int>(dMRI_image.imgDims[2] - 1)); // Adjust axial index (z-axis)
    update();
}

void MainGlWidget::handleCoronalScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120;
    j = std::clamp(j + delta, 0, static_cast<int>(dMRI_image.imgDims[1] - 1)); // Adjust coronal index (y-axis)
    update();
}

void MainGlWidget::handleSagittalScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120;
    i = std::clamp(i + delta, 0, static_cast<int>(dMRI_image.imgDims[0] - 1)); // Adjust sagittal index (x-axis)
    update();
}

void MainGlWidget::wheelEvent(QWheelEvent *event)
{
    // Check if Ctrl key is held
    if (event->modifiers() & Qt::ControlModifier) {
        // Handle zooming
        int delta = event->angleDelta().y();
        float zoomStep = 0.1f; // Adjust zoom step as needed
        if (delta > 0) {
            zoomFactor *= (1.0f + zoomStep);
        } else {
            zoomFactor *= (1.0f - zoomStep);
            if (zoomFactor < 0.1f) {
                zoomFactor = 0.1f; // Prevent zooming out too much
            }
        }
        update();
    } else {
        // Existing code to handle scrolling through slices
        QPointF mousePos = event->position();
        int widgetWidth = width();
        int viewportWidth = widgetWidth / 3;

        if (mousePos.x() < viewportWidth) {
            handleAxialScroll(event);
        } else if (mousePos.x() < 2 * viewportWidth) {
            handleCoronalScroll(event);
        } else {
            handleSagittalScroll(event);
        }
    }
}

void MainGlWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_R) {
        // Reset zoom and pan
        zoomFactor = 1.0f;
        panOffset = QVector3D(0.0f, 0.0f, 0.0f);
        update();
    }
}

void MainGlWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    update();  // Request a re-draw
}
