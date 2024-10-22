#include "mainglwidget.h"

MainGlWidget::MainGlWidget(QWidget *parent)
    : QOpenGLWidget(parent), dMRI_image("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/devices/MRI/T1.nii.gz"), mousePressed(false)
{
    QTimer *timer = new QTimer(this);

    std::cout << "Loading dMRI image" << std::endl;
    dMRI_image.read();
    dMRI_image.printInfo();
    std::cout << "Image voxel count: " << dMRI_image.voxCnt << std::endl;
    std::cout << "Number of dimensions:   " << dMRI_image.numberOfDimensions << std::endl;
    std::cout << "Dimensions:             [";
    for (int idx = 0; idx < dMRI_image.numberOfDimensions; idx++) {
        std::cout << dMRI_image.imgDims[idx];
        if (idx != (dMRI_image.numberOfDimensions - 1))
            std::cout << " x ";
    }
    std::cout << "]" << std::endl;

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
    auto setupOrtho = [](int imgWidth, int imgHeight) {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        // Standard coordinate system: origin at bottom-left
        glOrtho(0, imgWidth, 0, imgHeight, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    };

    // --------------------
    // Draw Axial Slice
    // --------------------
    glViewport(axialViewport.x, axialViewport.y, axialViewport.width, axialViewport.height);
    setupOrtho(dimX, dimY);

    int z = std::clamp(k, 0, dimZ - 1);

    // Draw the axial slice
    glBegin(GL_QUADS);
    for (int y = 0; y < dimY; ++y) {
        int flippedY = dimY - y - 1; // Flip the y-coordinate
        for (int x = 0; x < dimX; ++x) {
            int idx = voxelIndex(x, flippedY, dimZ - z - 1);
            float value = voxelData[idx];

            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            glColor3f(gray, gray, gray);

            float x0 = x;
            float y0 = y;
            float x1 = x + 1;
            float y1 = y + 1;

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
    glVertex2f(i, 0);
    glVertex2f(i, dimY);
    // Horizontal line at y = j
    float j_flipped = dimY - j - 1;
    glVertex2f(0, j_flipped);
    glVertex2f(dimX, j_flipped);
    glEnd();

    // --------------------
    // Draw Coronal Slice
    // --------------------
    glViewport(coronalViewport.x, coronalViewport.y, coronalViewport.width, coronalViewport.height);
    setupOrtho(dimX, dimZ);

    int y = std::clamp(j, 0, dimY - 1);

    // Draw the coronal slice
    glBegin(GL_QUADS);
    for (int z = 0; z < dimZ; ++z) {
        int flippedZ = dimZ - z - 1; // Flip the z-coordinate
        for (int x = 0; x < dimX; ++x) {
            int idx = voxelIndex(x, y, flippedZ);
            float value = voxelData[idx];

            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            glColor3f(gray, gray, gray);

            float x0 = x;
            float y0 = z;
            float x1 = x + 1;
            float y1 = z + 1;

            // Adjust y-coordinates to flip the image
            float y0_flipped = dimZ - y0 - 1;
            float y1_flipped = dimZ - y1 - 1;

            glVertex2f(x0, y0_flipped);
            glVertex2f(x1, y0_flipped);
            glVertex2f(x1, y1_flipped);
            glVertex2f(x0, y1_flipped);
        }
    }
    glEnd();

    // Draw lines on coronal slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    // Vertical line at x = i
    glVertex2f(i, 0);
    glVertex2f(i, dimZ);
    // Horizontal line at z = k (adjusted for flip)
    float k_flipped = dimZ - k - 1;
    glVertex2f(0, k_flipped);
    glVertex2f(dimX, k_flipped);
    glEnd();

    // --------------------
    // Draw Sagittal Slice
    // --------------------
    glViewport(sagittalViewport.x, sagittalViewport.y, sagittalViewport.width, widgetHeight);
    setupOrtho(dimY, dimZ);

    int x = std::clamp(i, 0, dimX - 1);

    // Draw the sagittal slice
    glBegin(GL_QUADS);
    for (int z = 0; z < dimZ; ++z) {
        int flippedZ = dimZ - z - 1; // Flip the z-coordinate
        for (int y = 0; y < dimY; ++y) {
            int idx = voxelIndex(x, y, flippedZ);
            float value = voxelData[idx];

            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            glColor3f(gray, gray, gray);

            float x0 = y;
            float y0 = z;
            float x1 = y + 1;
            float y1 = z + 1;

            // Adjust y-coordinates to flip the image
            float y0_flipped = dimZ - y0 - 1;
            float y1_flipped = dimZ - y1 - 1;

            glVertex2f(x0, y0_flipped);
            glVertex2f(x1, y0_flipped);
            glVertex2f(x1, y1_flipped);
            glVertex2f(x0, y1_flipped);
        }
    }
    glEnd();

    // Draw lines on sagittal slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    // Vertical line at y = j
    glVertex2f(j, 0);
    glVertex2f(j, dimZ);
    // Horizontal line at z = k (adjusted for flip)
    glVertex2f(0, k_flipped);
    glVertex2f(dimY, k_flipped);
    glEnd();
}
void MainGlWidget::handleAxialClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Calculate the position within the axial viewport
    int xInViewport = mousePos.x();
    int yInViewport = mousePos.y();

    // Map y-coordinate to OpenGL coordinate system (origin at bottom-left)
    int yFlipped = viewportHeight - yInViewport;

    // Dimensions of the Axial slice
    int dimX = dMRI_image.imgDims[0];
    int dimY = dMRI_image.imgDims[1];

    // Scale mouse coordinates to image coordinates
    float xRatio = static_cast<float>(dimX) / viewportWidth;
    float yRatio = static_cast<float>(dimY) / viewportHeight;

    int newI = static_cast<int>(xInViewport * xRatio);
    int newJ = static_cast<int>(yFlipped * yRatio);

    // Flip j due to image flipping
    newJ = dimY - newJ - 1;

    // Clamp the indices to valid ranges
    i = std::clamp(newI, 0, dimX - 1);
    j = std::clamp(newJ, 0, dimY - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::handleCoronalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Calculate the position within the coronal viewport
    int xInViewport = mousePos.x() - viewportWidth; // Subtract the offset
    int yInViewport = mousePos.y();

    int yFlipped = viewportHeight - yInViewport;

    // Dimensions of the Coronal slice
    int dimX = dMRI_image.imgDims[0];
    int dimZ = dMRI_image.imgDims[2];

    float xRatio = static_cast<float>(dimX) / viewportWidth;
    float zRatio = static_cast<float>(dimZ) / viewportHeight;

    int newI = static_cast<int>(xInViewport * xRatio);
    int newK = static_cast<int>(yFlipped * zRatio);

    // Flip k due to image flipping
    newK = dimZ - newK - 1;

    // Clamp the indices to valid ranges
    i = std::clamp(newI, 0, dimX - 1);
    k = std::clamp(newK, 0, dimZ - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::handleSagittalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Calculate the position within the sagittal viewport
    int xInViewport = mousePos.x() - 2 * viewportWidth; // Subtract the offset
    int yInViewport = mousePos.y();

    int yFlipped = viewportHeight - yInViewport;

    // Dimensions of the Sagittal slice
    int dimY = dMRI_image.imgDims[1];
    int dimZ = dMRI_image.imgDims[2];

    float yRatio = static_cast<float>(dimY) / viewportWidth;
    float zRatio = static_cast<float>(dimZ) / viewportHeight;

    int newJ = static_cast<int>(xInViewport * yRatio);
    int newK = static_cast<int>(yFlipped * zRatio);

    // Flip k due to image flipping
    newK = dimZ - newK - 1;

    // Clamp the indices to valid ranges
    j = std::clamp(newJ, 0, dimY - 1);
    k = std::clamp(newK, 0, dimZ - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::mousePressEvent(QMouseEvent *event)
{
    // Set the mouse pressed flag
    mousePressed = true;

    // Process the click as before
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

void MainGlWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (!mousePressed) {
        // If the mouse button is not pressed, do nothing
        return;
    }

    // Get the current mouse position
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
    // Get the position of the mouse cursor
    QPointF mousePos = event->position();

    int widgetWidth = width();
    int widgetHeight = height();
    int viewportWidth = widgetWidth / 3;

    // Determine which viewport the cursor is over
    if (mousePos.x() < viewportWidth) {
        // Cursor is over Axial Slice
        handleAxialScroll(event);
    } else if (mousePos.x() < 2 * viewportWidth) {
        // Cursor is over Coronal Slice
        handleCoronalScroll(event);
    } else {
        // Cursor is over Sagittal Slice
        handleSagittalScroll(event);
    }
}

void MainGlWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    update();  // Request a re-draw
}
