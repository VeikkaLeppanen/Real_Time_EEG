#include "mainglwidget.h"

MainGlWidget::MainGlWidget(QWidget *parent)
    : QOpenGLWidget(parent),
      mousePressed(false),
      zoomFactor(1.0f),
      panOffset(0.0f, 0.0f, 0.0f)
{
    // Initialize minValue and maxValue
    minValue = 0.0f;
    maxValue = 1.0f;

    // Initialize slice indices
    i = j = k = 0;

    overlayOpacity = 0.5f;
    
    this->setFocusPolicy ( Qt::StrongFocus );
    setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    
    // Set up the timer for updating the graph
    QTimer *timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &MainGlWidget::updateGraph);
    // timer->start(16); // Update approximately every 16 ms (60 FPS)
}

void MainGlWidget::setSliceIndices(int new_i, int new_j, int new_k)
{
    i = new_i;
    j = new_j;
    k = new_k;
    update(); // Trigger a repaint
}

Eigen::Matrix4f MainGlWidget::constructMatrix(float ijk2xyz[3][4]) {
    Eigen::Matrix4f mat = Eigen::Matrix4f::Identity();
    for (int i = 0; i < 3; ++i) {
        mat(i, 0) = ijk2xyz[i][0];
        mat(i, 1) = ijk2xyz[i][1];
        mat(i, 2) = ijk2xyz[i][2];
        mat(i, 3) = ijk2xyz[i][3];
    }
    return mat;
}

void MainGlWidget::loadImage_T1(const QString& filePath)
{
    // Load the new MRI image
    std::cout << "\nLoading MRI image from: " << filePath.toStdString() << std::endl;

    NIBR::Image<float> newImage(filePath.toStdString());
    T1_image = newImage;
    T1_image.read();
    T1_image.printInfo();

    // Check if the image was loaded successfully
    if (T1_image.data == nullptr || T1_image.voxCnt == 0) {
        std::cerr << "Failed to load image from: " << filePath.toStdString() << std::endl;
        return;
    }

    // Recalculate minValue and maxValue
    minValue = FLT_MAX;
    maxValue = -FLT_MAX;
    for (int idx = 0; idx < T1_image.voxCnt; ++idx) {
        float value = T1_image.data[idx];
        if (value < minValue) minValue = value;
        if (value > maxValue) maxValue = value;
    }
    std::cout << "Voxel value range: [" << minValue << ", " << maxValue << "]" << std::endl;

    // Update slice indices to the middle of the new image dimensions
    i = T1_image.imgDims[0] / 2; // Sagittal (x-axis)
    j = T1_image.imgDims[1] / 2; // Coronal (y-axis)
    k = T1_image.imgDims[2] / 2; // Axial (z-axis)

    T1_orientation = T1_image.getOrientation();
    std::cout << "Image orientation: " << T1_orientation[0] << T1_orientation[1] << T1_orientation[2] << std::endl;

    T1_pixDims[0] = T1_image.pixDims[0];
    T1_pixDims[1] = T1_image.pixDims[1];
    T1_pixDims[2] = T1_image.pixDims[2];

    // Trigger a repaint
    update();
}

void MainGlWidget::loadImage_fMRI(const QString& filePath)
{
    // Load the new MRI image
    std::cout << "\nLoading fMRI image from: " << filePath.toStdString() << std::endl;

    NIBR::Image<float> newImage(filePath.toStdString());
    fMRI_image = newImage;
    fMRI_image.read();
    fMRI_image.printInfo();

    // Check if the image was loaded successfully
    if (fMRI_image.data == nullptr || fMRI_image.voxCnt == 0) {
        std::cerr << "Failed to load image from: " << filePath.toStdString() << std::endl;
        return;
    }

    fMRI_image_resampled.createFromTemplate(T1_image, true);
    // imgResample(&fMRI_image_resampled, &fMRI_image);
    
    float* p = new float[3];
    for (int n = 0; n < fMRI_image_resampled.voxCnt; n++) {   
        fMRI_image_resampled.to_xyz(n,p);   
        fMRI_image_resampled.data[n] = fMRI_image(p, static_cast<int64_t>(0));
    }

    delete[] p;

    // Recalculate minValue and maxValue
    minValue_f = FLT_MAX;
    maxValue_f = -FLT_MAX;
    for (int idx = 0; idx < fMRI_image_resampled.voxCnt; ++idx) {
        float value = fMRI_image_resampled.data[idx];
        if (value < minValue_f) minValue_f = value;
        if (value > maxValue_f) maxValue_f = value;
    }
    
    std::cout << "Voxel value range: [" << minValue_f << ", " << maxValue_f << "]" << std::endl;

    // Trigger a repaint
    update();
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
    
    // Check if image data is available
    if (T1_image.voxCnt == 0) {
        // No image loaded; you might want to display a message or a blank screen
        return;
    }

    // Dimensions of the MRI image
    int dimX = T1_image.imgDims[0];
    int dimY = T1_image.imgDims[1];
    int dimZ = T1_image.imgDims[2];

    float pixDimX_T1 = T1_pixDims[0];//T1_image.pixDims[0];
    float pixDimY_T1 = T1_pixDims[1];//T1_image.pixDims[1];
    float pixDimZ_T1 = T1_pixDims[2];//T1_image.pixDims[2];

    int dimX_scaled = dimX * pixDimX_T1;
    int dimY_scaled = dimY * pixDimY_T1;
    int dimZ_scaled = dimZ * pixDimZ_T1;

    float* voxelData = T1_image.data;

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
    auto setupOrtho = [&](float imgWidth, float imgHeight, float panX, float panY, int viewportWidth, int viewportHeight) {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        float halfWidth = imgWidth / 2.0f;
        float halfHeight = imgHeight / 2.0f;

        // Calculate aspect ratios
        float aspectImage = imgWidth / imgHeight;
        float aspectViewport = static_cast<float>(viewportWidth) / viewportHeight;

        float left, right, bottom, top;

        if (aspectViewport > aspectImage) {
            // Viewport is wider than image; adjust width
            float adjustedHalfWidth = halfHeight * aspectViewport;
            left = (-adjustedHalfWidth - panX) / zoomFactor;
            right = (adjustedHalfWidth - panX) / zoomFactor;
            bottom = (-halfHeight - panY) / zoomFactor;
            top = (halfHeight - panY) / zoomFactor;
        } else {
            // Viewport is taller than image; adjust height
            float adjustedHalfHeight = halfWidth / aspectViewport;
            left = (-halfWidth - panX) / zoomFactor;
            right = (halfWidth - panX) / zoomFactor;
            bottom = (-adjustedHalfHeight - panY) / zoomFactor;
            top = (adjustedHalfHeight - panY) / zoomFactor;
        }

        glOrtho(left, right, bottom, top, -1, 1);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    };

    // --------------------
    // Draw Axial Slice
    // --------------------
    glViewport(axialViewport.x, axialViewport.y, axialViewport.width, axialViewport.height);

    // Calculate the physical dimensions of the image in the axial plane
    float imgWidth = dimX * pixDimX_T1;
    float imgHeight = dimY * pixDimY_T1;

    // Use the panOffset and zoomFactor
    setupOrtho(imgWidth, imgHeight, panOffset.x(), panOffset.y(), axialViewport.width, axialViewport.height);

    int z = std::clamp(k, 0, dimZ - 1);

    // Draw the T1 axial slice
    glBegin(GL_QUADS);
    for (int y = 0; y < dimY; ++y) {
        
        int yFlipped = y;
        if (T1_orientation[1] == "P") yFlipped = dimY - 1 - y;
        
        for (int x = 0; x < dimX; ++x) {
            
            int xFlipped = x;
            if (T1_orientation[0] == "R") xFlipped = dimX - 1 - x;
            
            int idx = voxelIndex(xFlipped, yFlipped, z);
            float value = voxelData[idx];

            // Normalize T1 voxel value for grayscale color
            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            // Base color from T1 image
            float r = gray;
            float g = gray;
            float b = gray;
            float a = 1.0f; // Full opacity for T1 image

            // Initialize overlay intensity
            float overlayIntensity = 0.0f;

            // Check if fMRI data is available
            if (fMRI_image_resampled.voxCnt != 0) {
                // If images are resampled and aligned, use the same index
                float value_f = fMRI_image_resampled.data[idx];

                // Normalize fMRI value for visualization
                float intensity = std::clamp((value_f - minValue_f) / (maxValue_f - minValue_f), 0.0f, 1.0f);

                // Update overlay intensity based on the fMRI value
                overlayIntensity = intensity * overlayOpacity;
            }

            // Combine the T1 and fMRI colors
            // For example, enhance the red channel based on the overlay
            r = r * (1.0f - overlayIntensity) + overlayIntensity;

            // Set the final color with alpha blending
            glColor4f(r, g, b, a);

            // Adjust coordinates to center around (0,0)
            float x0 = (x * pixDimX_T1) - (imgWidth / 2.0f);
            float y0 = (y * pixDimY_T1) - (imgHeight / 2.0f);
            float x1 = x0 + pixDimX_T1;
            float y1 = y0 + pixDimY_T1;

            // Draw the quad
            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();

    // Disable blending after overlay
    glDisable(GL_BLEND);

    // Draw crosshairs on axial slice
    // Draw crosshairs on axial slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);

    // Adjust x index for orientation
    int xCrosshair = i;
    if (T1_orientation[0] == "R") xCrosshair = dimX - 1 - i;

    // Compute lineX in physical coordinates
    float lineX = (xCrosshair * pixDimX_T1) - (imgWidth / 2.0f) + (pixDimX_T1 / 2.0f);

    // Adjust y index for orientation
    int yCrosshair = j;
    if (T1_orientation[1] == "P") yCrosshair = dimY - 1 - j;

    // Compute lineY in physical coordinates
    float lineY = (yCrosshair * pixDimY_T1) - (imgHeight / 2.0f) + (pixDimY_T1 / 2.0f);

    // Vertical line at x = lineX
    glVertex2f(lineX, -imgHeight / 2.0f);
    glVertex2f(lineX, imgHeight / 2.0f);

    // Horizontal line at y = lineY
    glVertex2f(-imgWidth / 2.0f, lineY);
    glVertex2f(imgWidth / 2.0f, lineY);

    glEnd();

    // --------------------
    // Draw Coronal Slice
    // --------------------
    glViewport(coronalViewport.x, coronalViewport.y, coronalViewport.width, widgetHeight);

    // Calculate the physical dimensions of the image in the coronal plane
    imgWidth = dimX * pixDimX_T1;
    imgHeight = dimZ * pixDimZ_T1;

    // Use the panOffset and zoomFactor
    setupOrtho(imgWidth, imgHeight, panOffset.x(), panOffset.z(), coronalViewport.width, coronalViewport.height);

    int y = std::clamp(j, 0, dimY - 1);

    // Draw the T1 coronal slice
    glBegin(GL_QUADS);
    for (int z = 0; z < dimZ; ++z) {

        int zFlipped = z;
        if (T1_orientation[2] == "I") zFlipped = dimZ - 1 - z;

        for (int x = 0; x < dimX; ++x) {

            int xFlipped = x;
            if (T1_orientation[0] == "R") xFlipped = dimX - 1 - x;

            int idx = voxelIndex(xFlipped, y, zFlipped);
            float value = voxelData[idx];

            // Normalize T1 voxel value for grayscale color
            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            // Base color from T1 image
            float r = gray;
            float g = gray;
            float b = gray;
            float a = 1.0f; // Full opacity for T1 image

            // Initialize overlay intensity
            float overlayIntensity = 0.0f;

            // Check if fMRI data is available
            if (fMRI_image_resampled.voxCnt != 0) {
                // If images are resampled and aligned, use the same index
                float value_f = fMRI_image_resampled.data[idx];

                // Normalize fMRI value for visualization
                float intensity = std::clamp((value_f - minValue_f) / (maxValue_f - minValue_f), 0.0f, 1.0f);

                // Update overlay intensity based on the fMRI value
                overlayIntensity = intensity * overlayOpacity;
            }

            // Combine the T1 and fMRI colors
            // For example, enhance the red channel based on the overlay
            r = r * (1.0f - overlayIntensity) + overlayIntensity;

            // Set the final color with alpha blending
            glColor4f(r, g, b, a);

            // Adjust coordinates to center around (0,0)
            float x0 = (x * pixDimX_T1) - (imgWidth / 2.0f);
            float y0 = (z * pixDimZ_T1) - (imgHeight / 2.0f);
            float x1 = x0 + pixDimX_T1;
            float y1 = y0 + pixDimZ_T1;

            // Draw the quad
            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();

    // Draw crosshairs on axial slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    
    // Adjust x index for orientation
    xCrosshair = i;
    if (T1_orientation[0] == "R") xCrosshair = dimX - 1 - i;
    
    // Compute lineX in physical coordinates
    lineX = (xCrosshair * pixDimX_T1) - (imgWidth / 2.0f) + (pixDimX_T1 / 2.0f);
    
    // Adjust z index for orientation
    yCrosshair = k;
    if (T1_orientation[2] == "I") yCrosshair = dimZ - 1 - k;
    
    // Compute lineY in physical coordinates
    lineY = (yCrosshair * pixDimZ_T1) - (imgHeight / 2.0f) + (pixDimZ_T1 / 2.0f);
    
    // Vertical line at x = lineX
    glVertex2f(lineX, -imgHeight / 2.0f);
    glVertex2f(lineX, imgHeight / 2.0f);
    
    // Horizontal line at z = lineY
    glVertex2f(-imgWidth / 2.0f, lineY);
    glVertex2f(imgWidth / 2.0f, lineY);
    
    glEnd();

    // --------------------
    // Draw Sagittal Slice
    // --------------------
    glViewport(sagittalViewport.x, sagittalViewport.y, sagittalViewport.width, widgetHeight);

    // Calculate the physical dimensions of the image in the sagittal plane
    imgWidth = dimY * pixDimY_T1;
    imgHeight = dimZ * pixDimZ_T1;
    
    // Use the panOffset and zoomFactor
    setupOrtho(imgWidth, imgHeight, -panOffset.y(), panOffset.z(), sagittalViewport.width, sagittalViewport.height);

    int x = std::clamp(i, 0, dimX - 1);

    // Draw the T1 sagittal slice
    glBegin(GL_QUADS);
    for (int z = 0; z < dimZ; ++z) {

        int zFlipped = z;
        if (T1_orientation[2] == "I") zFlipped = dimZ - 1 - z;

        for (int y = 0; y < dimY; ++y) {

            int yFlipped = y;
            if (T1_orientation[1] == "A") yFlipped = dimY - 1 - y;

            int idx = voxelIndex(x, yFlipped, zFlipped);
            float value = voxelData[idx];

            // Normalize T1 voxel value for grayscale color
            float gray = std::clamp((value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

            // Base color from T1 image
            float r = gray;
            float g = gray;
            float b = gray;
            float a = 1.0f; // Full opacity for T1 image

            // Initialize overlay intensity
            float overlayIntensity = 0.0f;

            // Check if fMRI data is available
            if (fMRI_image_resampled.voxCnt != 0) {
                // If images are resampled and aligned, use the same index
                float value_f = fMRI_image_resampled.data[idx];

                // Normalize fMRI value for visualization
                float intensity = std::clamp((value_f - minValue_f) / (maxValue_f - minValue_f), 0.0f, 1.0f);

                // Update overlay intensity based on the fMRI value
                overlayIntensity = intensity * overlayOpacity;
            }

            // Combine the T1 and fMRI colors
            // For example, enhance the red channel based on the overlay
            r = r * (1.0f - overlayIntensity) + overlayIntensity;

            // Set the final color with alpha blending
            glColor4f(r, g, b, a);

            // Adjust coordinates to center around (0,0)
            float x0 = (y * pixDimY_T1) - (imgWidth / 2.0f);
            float y0 = (z * pixDimZ_T1) - (imgHeight / 2.0f);
            float x1 = x0 + pixDimY_T1;
            float y1 = y0 + pixDimZ_T1;

            // Draw the quad
            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
        }
    }
    glEnd();

    // Draw crosshairs on sagittal slice
    glColor3f(0.5f, 0.5f, 0.0f); // Yellow color
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    
    // Adjust x index for orientation
    xCrosshair = j;
    if (T1_orientation[1] == "A") xCrosshair = dimY - 1 - j;
    
    // Compute lineX in physical coordinates
    lineX = (xCrosshair * pixDimY_T1) - (imgWidth / 2.0f) + (pixDimY_T1 / 2.0f);
    
    // Adjust z index for orientation
    yCrosshair = k;
    if (T1_orientation[2] == "I") yCrosshair = dimZ - 1 - k;
    
    // Compute lineY in physical coordinates
    lineY = (yCrosshair * pixDimZ_T1) - (imgHeight / 2.0f) + (pixDimZ_T1 / 2.0f);
    
    // Vertical line at x = lineX
    glVertex2f(lineX, -imgHeight / 2.0f);
    glVertex2f(lineX, imgHeight / 2.0f);
    
    // Horizontal line at z = lineY
    glVertex2f(-imgWidth / 2.0f, lineY);
    glVertex2f(imgWidth / 2.0f, lineY);
    
    glEnd();
    
    // Retrieve the voxel indices
    int i_trgt = std::clamp(i, 0, dimX - 1);
    int j_trgt = std::clamp(j, 0, dimY - 1);
    int k_trgt = std::clamp(k, 0, dimZ - 1);


    // Retrieve the target voxel value
    int trgt_idx = voxelIndex(i_trgt, j_trgt, k_trgt);
    float trgt_voxelValue = voxelData[trgt_idx];
    float trgt_voxelValue_f = 0;

    // Retrieve the voxel positions
    float p[3];
    T1_image.to_xyz(trgt_idx, p);

    if (fMRI_image.voxCnt != 0) {
        // Retrieve the target voxel value
        trgt_voxelValue_f = fMRI_image_resampled.data[trgt_idx];
    }

    // Begin QPainter
    QPainter painter(this);

    // Set text properties
    painter.setPen(Qt::yellow);
    painter.setFont(QFont("Arial", 10));

    // Get the position unit
    QString spaceUnit = QString::fromStdString(T1_image.getSpaceUnit());

    // Prepare the text
    QString text = QString("value T1: %1\nvalue fMRI: %2\nvoxel index: [ %3, %4, %5 ]\nposition: [ %6, %7, %8 ] %9")
        .arg(trgt_voxelValue)
        .arg(trgt_voxelValue_f)
        .arg(i_trgt).arg(j_trgt).arg(k_trgt)
        .arg(p[0]).arg(p[1]).arg(p[2]).arg(spaceUnit);

    // Position the text at the bottom left corner
    QRect textRect(10, height() - 110, 300, 100); // x, y, width, height

    painter.drawText(textRect, Qt::AlignLeft | Qt::AlignBottom | Qt::TextWordWrap, text);

    // End QPainter
    painter.end();
}

void MainGlWidget::handleAxialClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Dimensions of the Axial slice
    int dimX = T1_image.imgDims[0];
    int dimY = T1_image.imgDims[1];

    // Calculate the position within the axial viewport
    int xInViewport = viewportWidth - mousePos.x();
    int yInViewport = viewportHeight - mousePos.y(); // Flip y-coordinate

    // Convert viewport coordinates to Normalized Device Coordinates (NDC)
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewport) / viewportHeight) * 2.0f - 1.0f;

    // Undo the projection transformation (glOrtho)
    float halfWidth = dimX / 2.0f;
    float halfHeight = dimY / 2.0f;

    float panX = panOffset.x();
    if (T1_orientation[0] == "R") panX = -panOffset.x();
    panX = panX / T1_pixDims[0];
    float panY = panOffset.y();
    if (T1_orientation[1] == "P") panY = -panOffset.y();
    panY = panY / T1_pixDims[1];

    float left = (-halfWidth - panX) / zoomFactor;
    float right = (halfWidth - panX) / zoomFactor;
    float bottom = (-halfHeight - panY) / zoomFactor;
    float top = (halfHeight - panY) / zoomFactor;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (right - left) + left;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (top - bottom) + bottom;

    // Convert world coordinates to image coordinates
    int imgX = static_cast<int>(worldX + dimX / 2.0f);
    int imgY = static_cast<int>(worldY + dimY / 2.0f);

    // Clamp the indices to valid ranges
    i = std::clamp(imgX, 0, dimX - 1);
    j = std::clamp(imgY, 0, dimY - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::handleCoronalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Dimensions of the Coronal slice
    int dimX = T1_image.imgDims[0];
    int dimZ = T1_image.imgDims[2];

    // Adjust mouse position for the coronal viewport
    int xInViewport = (viewportWidth - mousePos.x()) + viewportWidth; // Subtract viewport offset
    int yInViewport = viewportHeight - mousePos.y(); // Flip y-coordinate

    // Convert viewport coordinates to NDC
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewport) / viewportHeight) * 2.0f - 1.0f;

    // Undo the projection transformation
    float halfWidth = dimX / 2.0f;
    float halfHeight = dimZ / 2.0f;

    float panX = panOffset.x();
    if (T1_orientation[0] == "R") panX = -panOffset.x();
    panX = panX / T1_pixDims[0];
    float panZ = panOffset.z();
    if (T1_orientation[2] == "I") panZ = -panOffset.z();
    panZ = panZ / T1_pixDims[2];

    float left = (-halfWidth - panX) / zoomFactor;
    float right = (halfWidth - panX) / zoomFactor;
    float bottom = (-halfHeight - panZ) / zoomFactor;
    float top = (halfHeight - panZ) / zoomFactor;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (right - left) + left;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (top - bottom) + bottom;

    // Convert world coordinates to image coordinates
    int imgX = static_cast<int>(worldX + dimX / 2.0f);
    int imgZ = static_cast<int>(worldY + dimZ / 2.0f);

    // Clamp the indices to valid ranges
    i = std::clamp(imgX, 0, dimX - 1);
    k = std::clamp(imgZ, 0, dimZ - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::handleSagittalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    // Dimensions of the Sagittal slice
    int dimY = T1_image.imgDims[1];
    int dimZ = T1_image.imgDims[2];

    // Adjust mouse position for the sagittal viewport
    int xInViewport = (viewportWidth - mousePos.x()) + 2 * viewportWidth; // Subtract viewport offset
    int yInViewport = viewportHeight - mousePos.y(); // Flip y-coordinate

    // Convert viewport coordinates to NDC
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewport) / viewportHeight) * 2.0f - 1.0f;

    // Undo the projection transformation
    float halfWidth = dimY / 2.0f;
    float halfHeight = dimZ / 2.0f;

    float panY = -panOffset.y();
    if (T1_orientation[1] == "A") panY = panOffset.y();
    panY = panY / T1_pixDims[1];
    float panZ = panOffset.z();
    if (T1_orientation[2] == "I") panZ = -panOffset.z();
    panZ = panZ / T1_pixDims[2];

    float left = (-halfWidth - panY) / zoomFactor;
    float right = (halfWidth - panY) / zoomFactor;
    float bottom = (-halfHeight - panZ) / zoomFactor;
    float top = (halfHeight - panZ) / zoomFactor;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (right - left) + left;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (top - bottom) + bottom;

    // Convert world coordinates to image coordinates
    int imgY = static_cast<int>(worldX + dimY / 2.0f);
    int imgZ = static_cast<int>(worldY + dimZ / 2.0f);

    // Clamp the indices to valid ranges
    j = std::clamp(imgY, 0, dimY - 1);
    k = std::clamp(imgZ, 0, dimZ - 1);

    update(); // Trigger a repaint
}

void MainGlWidget::mousePressEvent(QMouseEvent *event)
{
    if (T1_image.data == nullptr || T1_image.voxCnt == 0) {
        // No image loaded; ignore the event
        return;
    }

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
    if (T1_image.data == nullptr || T1_image.voxCnt == 0) {
        // No image loaded; ignore the event
        return;
    }

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
            scaleX = (T1_image.imgDims[0] * T1_pixDims[0] + zoomFactor) / viewportWidth;
            scaleY = (T1_image.imgDims[1] * T1_pixDims[1] + zoomFactor) / widgetHeight;

            panOffset.setX((panOffset.x() - deltaX * scaleX));
            panOffset.setY((panOffset.y() - deltaY * scaleY));

        } else if (viewportIndex == 1) {
            // Coronal slice (x, z plane)
            scaleX = (T1_image.imgDims[0] * T1_pixDims[0] + zoomFactor) / viewportWidth;
            scaleY = (T1_image.imgDims[2] * T1_pixDims[2] + zoomFactor) / widgetHeight;

            panOffset.setX(panOffset.x() - deltaX * scaleX);
            panOffset.setZ(panOffset.z() - deltaY * scaleY);

        } else {
            // Sagittal slice (y, z plane)
            scaleX = (T1_image.imgDims[1] * T1_pixDims[1] + zoomFactor) / viewportWidth;
            scaleY = (T1_image.imgDims[2] * T1_pixDims[2] + zoomFactor) / widgetHeight;

            panOffset.setY(panOffset.y() + deltaX * scaleX);
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
    if (T1_image.data == nullptr || T1_image.voxCnt == 0) {
        // No image loaded; ignore the event
        return;
    }

    Q_UNUSED(event);
    // Reset the mouse pressed flag
    mousePressed = false;
}

void MainGlWidget::handleAxialScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120; // Each notch is 15 degrees; 120 units per notch
    k = std::clamp(k + delta, 0, static_cast<int>(T1_image.imgDims[2] - 1)); // Adjust axial index (z-axis)
    update();
}

void MainGlWidget::handleCoronalScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120;
    j = std::clamp(j + delta, 0, static_cast<int>(T1_image.imgDims[1] - 1)); // Adjust coronal index (y-axis)
    update();
}

void MainGlWidget::handleSagittalScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120;
    i = std::clamp(i + delta, 0, static_cast<int>(T1_image.imgDims[0] - 1)); // Adjust sagittal index (x-axis)
    update();
}

void MainGlWidget::wheelEvent(QWheelEvent *event)
{
    if (T1_image.data == nullptr || T1_image.voxCnt == 0) {
        // No image loaded; ignore the event
        return;
    }

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

// Add this function to perform trilinear interpolation
float MainGlWidget::getInterpolatedVoxelValue(float* data, float x, float y, float z, int dimX, int dimY, int dimZ) {
    int x0 = static_cast<int>(floor(x));
    int y0 = static_cast<int>(floor(y));
    int z0 = static_cast<int>(floor(z));

    int x1 = x0 + 1;
    int y1 = y0 + 1;
    int z1 = z0 + 1;

    float xd = x - x0;
    float yd = y - y0;
    float zd = z - z0;

    // Clamp indices to image dimensions
    x0 = std::clamp(x0, 0, dimX - 1);
    x1 = std::clamp(x1, 0, dimX - 1);
    y0 = std::clamp(y0, 0, dimY - 1);
    y1 = std::clamp(y1, 0, dimY - 1);
    z0 = std::clamp(z0, 0, dimZ - 1);
    z1 = std::clamp(z1, 0, dimZ - 1);

    // Calculate voxel indices directly
    int idx000 = x0 + y0 * dimX + z0 * dimX * dimY;
    int idx100 = x1 + y0 * dimX + z0 * dimX * dimY;
    int idx010 = x0 + y1 * dimX + z0 * dimX * dimY;
    int idx110 = x1 + y1 * dimX + z0 * dimX * dimY;
    int idx001 = x0 + y0 * dimX + z1 * dimX * dimY;
    int idx101 = x1 + y0 * dimX + z1 * dimX * dimY;
    int idx011 = x0 + y1 * dimX + z1 * dimX * dimY;
    int idx111 = x1 + y1 * dimX + z1 * dimX * dimY;

    // Retrieve voxel values at the eight corners
    float c000 = data[idx000];
    float c100 = data[idx100];
    float c010 = data[idx010];
    float c110 = data[idx110];
    float c001 = data[idx001];
    float c101 = data[idx101];
    float c011 = data[idx011];
    float c111 = data[idx111];

    // Perform trilinear interpolation
    float c00 = c000 * (1 - xd) + c100 * xd;
    float c01 = c001 * (1 - xd) + c101 * xd;
    float c10 = c010 * (1 - xd) + c110 * xd;
    float c11 = c011 * (1 - xd) + c111 * xd;

    float c0 = c00 * (1 - yd) + c10 * yd;
    float c1 = c01 * (1 - yd) + c11 * yd;

    float c = c0 * (1 - zd) + c1 * zd;

    return c;
}

void MainGlWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    update();  // Request a re-draw
}













// ----------------------------
// ROI functions
// ----------------------------
void MainGlWidget::addButton_clicked() 
{
    NIBR::Image<bool> newImage;
    newImage.createFromTemplate(T1_image, true);

    // Find the next available number for "ROI"
    int maxNumber = 0;
    for (const auto& name : ROI_names) {
        if (name.rfind("ROI", 0) == 0) { // Check if the name starts with "ROI"
            try {
                int number = std::stoi(name.substr(3)); // Extract number after "ROI"
                maxNumber = std::max(maxNumber, number);
            } catch (...) {
                // Ignore any non-numeric entries
            }
        }
    }

    std::string newName = "ROI" + std::to_string(maxNumber + 1);
    ROI_names.push_back(newName); // Add the new name
    ROI_vector.push_back(newImage);
    ROI_visibility.push_back(true);

    emit ROI_update(ROI_names, ROI_visibility);
}


void MainGlWidget::loadButton_clicked(const QString& filePath) 
{
    std::string filePath_std = filePath.toStdString();
    NIBR::Image<bool> newImage(filePath_std);
    newImage.read();
    ROI_vector.push_back(newImage);

    // Find the file name in filePath
    size_t lastSlash = filePath_std.find_last_of("/\\");
    size_t lastDot = filePath_std.find_last_of(".");
    std::string fileName = filePath_std.substr(lastSlash + 1, lastDot - lastSlash - 1);

    ROI_names.push_back(fileName);
    ROI_visibility.push_back(true);

    emit ROI_update(ROI_names, ROI_visibility);
}

void MainGlWidget::saveButton_clicked(const std::vector<bool>& states) 
{
    for (int i = 0; i < states.size(); i++) {
        if (states[i]) {
            ROI_vector[i].write(ROI_names[i] + ".nii.gz");
        }
    }
}

void MainGlWidget::deleteButton_clicked(const std::vector<bool>& states) 
{
    for (int i = 0; i < states.size(); i++) {
        if (states[i]) {
            ROI_names.erase(ROI_names.begin() + i);
            ROI_vector.erase(ROI_vector.begin() + i);
            ROI_visibility.erase(ROI_visibility.begin() + i);
        }
    }
    emit ROI_update(ROI_names, ROI_visibility);
}

void MainGlWidget::visibleButton_clicked(const std::vector<bool>& states) 
{
    for (int i = 0; i < states.size(); i++) {
        if (states[i]) {
            ROI_visibility[i] = !ROI_visibility[i];
        }
    }
    emit ROI_update(ROI_names, ROI_visibility);
}
