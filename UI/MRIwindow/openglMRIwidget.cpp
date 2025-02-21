#include "openglMRIwidget.h"

OpenGLMRIWidget::OpenGLMRIWidget(QWidget *parent)
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

    fMRI_overlayOpacity = 0.5f;
    
    this->setFocusPolicy ( Qt::StrongFocus );
    setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    
    // Set up the timer for updating the graph
    QTimer *timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, &OpenGLMRIWidget::updateGraph);
    // timer->start(16); // Update approximately every 16 ms (60 FPS)

    // Initialize file watcher
    fileWatcher = new QFileSystemWatcher(this);
    connect(fileWatcher, &QFileSystemWatcher::directoryChanged,
            this, &OpenGLMRIWidget::handleNewFile);
    setWatchFolder("/home/veikka/Work/EEG/DataStream/Real_Time_EEG/devices/MRI/Images");
}

void OpenGLMRIWidget::setWatchFolder(const QString& path)
{
    // Remove any existing watched paths
    if (!fileWatcher->directories().isEmpty()) {
        fileWatcher->removePaths(fileWatcher->directories());
    }
    if (!fileWatcher->files().isEmpty()) {
        fileWatcher->removePaths(fileWatcher->files());
    }
    
    // Set new watch path
    watchFolderPath = path;
    fileWatcher->addPath(path);
    std::cout << "Now watching folder: " << path.toStdString() << std::endl;

    // Store initial list of files
    QDir dir(path);
    QStringList filters;
    filters << "*.nii" << "*.nii.gz" << "*.dcm" << "*.dicom" << "*.ima" << "*.IMA" << "*";
    dir.setNameFilters(filters);
    dir.setFilter(QDir::Files | QDir::NoDotAndDotDot);
    processedFiles = dir.entryList();
}

void OpenGLMRIWidget::handleNewFile(const QString& path)
{
    // Get list of files in directory
    QDir dir(path);
    QStringList filters;
    filters << "*.nii" << "*.nii.gz" << "*.dcm" << "*.dicom" << "*.ima" << "*.IMA" << "*";
    dir.setNameFilters(filters);
    dir.setFilter(QDir::Files | QDir::NoDotAndDotDot);
    QStringList currentFiles = dir.entryList();
    
    // Find new files by comparing with processed files
    QStringList newFiles;
    for (const QString& file : currentFiles) {
        if (!processedFiles.contains(file)) {
            newFiles.append(file);
        }
    }

    // If there are new files, process the most recent one
    if (!newFiles.isEmpty()) {
        // Find the most recent file
        QString mostRecentFile;
        QDateTime mostRecentTime;

        for (const QString& file : newFiles) {
            QString fullPath = dir.filePath(file);
            QFileInfo fileInfo(fullPath);
            QDateTime created = fileInfo.birthTime();
            if (created.isValid() && (mostRecentFile.isEmpty() || created > mostRecentTime)) {
                mostRecentFile = file;
                mostRecentTime = created;
            }
        }

        if (!mostRecentFile.isEmpty()) {
            QString fullPath = dir.filePath(mostRecentFile);
            QString extension = QFileInfo(mostRecentFile).suffix().toLower();
            
            // Load file if it's a NIFTI, DICOM, or has no extension
            if (extension.isEmpty() || 
                extension == "nii" || 
                extension == "gz" || 
                extension == "dcm" || 
                extension == "dicom" || 
                extension == "ima") {
                std::cout << "Loading new fMRI file: " << fullPath.toStdString() << std::endl;
                loadImage_fMRI(fullPath);
            }
        }
    }

    // Update the list of processed files
    processedFiles = currentFiles;

    // Watch all files in the directory
    for (const QString& file : currentFiles) {
        QString fullPath = dir.filePath(file);
        if (!fileWatcher->files().contains(fullPath)) {
            fileWatcher->addPath(fullPath);
        }
    }
}

void OpenGLMRIWidget::setSliceIndices(int new_i, int new_j, int new_k)
{
    i = new_i;
    j = new_j;
    k = new_k;
    update(); // Trigger a repaint
}

Eigen::Matrix4f OpenGLMRIWidget::constructMatrix(float ijk2xyz[3][4]) {
    Eigen::Matrix4f mat = Eigen::Matrix4f::Identity();
    for (int i = 0; i < 3; ++i) {
        mat(i, 0) = ijk2xyz[i][0];
        mat(i, 1) = ijk2xyz[i][1];
        mat(i, 2) = ijk2xyz[i][2];
        mat(i, 3) = ijk2xyz[i][3];
    }
    return mat;
}

void OpenGLMRIWidget::loadImage_T1(const QString& filePath)
{
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
    mean_intensity = 0;

    for (int idx = 0; idx < T1_image.voxCnt; ++idx) {
        float value = T1_image.data[idx];
        if (value < minValue) minValue = value;
        if (value > maxValue) maxValue = value;
        mean_intensity += value;
    }
    mean_intensity /= T1_image.voxCnt;

    // Update slice indices to the middle of the new image dimensions
    i = T1_image.imgDims[0] / 2;
    j = T1_image.imgDims[1] / 2;
    k = T1_image.imgDims[2] / 2;

    T1_orientation = T1_image.getOrientation();
    T1_pixDims[0] = T1_image.pixDims[0];
    T1_pixDims[1] = T1_image.pixDims[1];
    T1_pixDims[2] = T1_image.pixDims[2];
}

void OpenGLMRIWidget::loadImage_fMRI(const QString& filePath)
{
    std::cout << "\nLoading fMRI image from: " << filePath.toStdString() << std::endl;

    NIBR::Image<float> newImage(filePath.toStdString());

    fMRI_image = newImage;
    fMRI_image.read();

    if (fMRI_image.data == nullptr || fMRI_image.voxCnt == 0) {
        std::cerr << "Failed to load image from: " << filePath.toStdString() << std::endl;
        return;
    }

            fMRI_image_resampled.createFromTemplate(T1_image, true);
            fMRI_image_resampled.createFromTemplate(T1_image, true);

    fMRI_image_resampled.createFromTemplate(T1_image, true);

    float* p = new float[3];

    for (int n = 0; n < fMRI_image_resampled.voxCnt; n++) {   
        fMRI_image_resampled.to_xyz(n,p);   
        fMRI_image_resampled.data[n] = fMRI_image(p, static_cast<int64_t>(0));
    }

    delete[] p;

    minValue_f = FLT_MAX;
    maxValue_f = -FLT_MAX;
    mean_intensity_f = 0;

    for (int idx = 0; idx < fMRI_image_resampled.voxCnt; ++idx) {
        float value = fMRI_image_resampled.data[idx];
        if (value < minValue_f) minValue_f = value;
        if (value > maxValue_f) maxValue_f = value;
        mean_intensity_f += value;
    }

    mean_intensity_f /= fMRI_image_resampled.voxCnt;

    calculateROIMeans();
}

void OpenGLMRIWidget::initializeGL()
{
    initializeOpenGLFunctions();
    glClearColor(0.0, 0.0, 0.0, 1.0); // Set the clearing color to black
}

void OpenGLMRIWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h); // Set the viewport to cover the entire widget
}

void OpenGLMRIWidget::paintGL()
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

    float pixDimX_T1 = T1_pixDims[0];
    float pixDimY_T1 = T1_pixDims[1];
    float pixDimZ_T1 = T1_pixDims[2];

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

    Viewport axialViewport = {0, 0, viewportWidth, widgetHeight};
    Viewport coronalViewport = {viewportWidth, 0, viewportWidth, widgetHeight};
    Viewport sagittalViewport = {2 * viewportWidth, 0, viewportWidth, widgetHeight};

    // Helper function to set up projection and modelview matrices
    auto setupOrtho = [&](float imgWidth, float imgHeight, float panX, float panY,
                          int viewportWidth, int viewportHeight,
                          float& left, float& right, float& bottom, float& top) {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        float halfWidth = imgWidth / 2.0f;
        float halfHeight = imgHeight / 2.0f;

        // Calculate aspect ratios
        float aspectImage = imgWidth / imgHeight;
        float aspectViewport = static_cast<float>(viewportWidth) / viewportHeight;

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
    imgWidth_axial = dimX * pixDimX_T1;
    imgHeight_axial = dimY * pixDimY_T1;

    // Use the panOffset and zoomFactor
    // For Axial Slice
    setupOrtho(imgWidth_axial, imgHeight_axial, panOffset.x(), panOffset.y(),
               axialViewport.width, axialViewport.height,
               leftAxial, rightAxial, bottomAxial, topAxial);

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

            // Adjust the voxel value using contrast factor
            float adjusted_value = (value - mean_intensity) * T1_contrast + mean_intensity;

            // Normalize T1 voxel value for grayscale color
            float gray = std::clamp((adjusted_value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

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

                // Adjust the voxel value using contrast factor
                float adjusted_value_f = (value_f - mean_intensity_f) * fMRI_contrast + mean_intensity_f;

                // Normalize fMRI value for visualization
                float intensity = std::clamp((adjusted_value_f - minValue_f) / (maxValue_f - minValue_f), 0.0f, 1.0f);

                // Update overlay intensity based on the fMRI value
                overlayIntensity = intensity * fMRI_overlayOpacity;
            }

            // Combine the T1 and fMRI colors
            r = r * (1.0f - overlayIntensity) + overlayIntensity;

            for (int ROI_num = 0; ROI_num < ROI_vector.size(); ROI_num++) {
                if(!ROI_toggle_states[ROI_num]) continue;

                NIBR::Image<bool> &ROI_image = *(ROI_vector[ROI_num]);
                QColor ROI_color = ROI_colors[ROI_num];
                float ROI_opacity = ROI_opacities[ROI_num];
                if (ROI_image.data[idx]) {
                    // Extract RGB components from the ROI color
                    float roi_r = ROI_color.redF();   // Normalize to range [0, 1]
                    float roi_g = ROI_color.greenF(); // Normalize to range [0, 1]
                    float roi_b = ROI_color.blueF();  // Normalize to range [0, 1]

                    // Blend the current color with the ROI color using alpha blending (ROI_opacity)
                    r = r * (1.0f - ROI_opacity) + roi_r * ROI_opacity;
                    g = g * (1.0f - ROI_opacity) + roi_g * ROI_opacity;
                    b = b * (1.0f - ROI_opacity) + roi_b * ROI_opacity;

                    // Adjust the alpha value (transparency) based on ROI_opacity
                    a = a * (1.0f - ROI_opacity) + ROI_opacity;
                }
            }

            // Set the final color with alpha blending
            glColor4f(r, g, b, a);

            // Adjust coordinates to center around (0,0)
            float x0 = (x * pixDimX_T1) - (imgWidth_axial / 2.0f);
            float y0 = (y * pixDimY_T1) - (imgHeight_axial / 2.0f);
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
    float lineX = (xCrosshair * pixDimX_T1) - (imgWidth_axial / 2.0f) + (pixDimX_T1 / 2.0f);

    // Adjust y index for orientation
    int yCrosshair = j;
    if (T1_orientation[1] == "P") yCrosshair = dimY - 1 - j;

    // Compute lineY in physical coordinates
    float lineY = (yCrosshair * pixDimY_T1) - (imgHeight_axial / 2.0f) + (pixDimY_T1 / 2.0f);

    // Vertical line at x = lineX
    glVertex2f(lineX, -imgHeight_axial / 2.0f);
    glVertex2f(lineX, imgHeight_axial / 2.0f);

    // Horizontal line at y = lineY
    glVertex2f(-imgWidth_axial / 2.0f, lineY);
    glVertex2f(imgWidth_axial / 2.0f, lineY);

    glEnd();

    // --------------------
    // Draw Coronal Slice
    // --------------------
    glViewport(coronalViewport.x, coronalViewport.y, coronalViewport.width, widgetHeight);

    // Calculate the physical dimensions of the image in the coronal plane
    imgWidth_coronal = dimX * pixDimX_T1;
    imgHeight_coronal = dimZ * pixDimZ_T1;

    // Use the panOffset and zoomFactor
    setupOrtho(imgWidth_coronal, imgHeight_coronal, panOffset.x(), panOffset.z(),
               coronalViewport.width, coronalViewport.height,
               leftCoronal, rightCoronal, bottomCoronal, topCoronal);
    
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

            // Adjust the voxel value using contrast factor
            float adjusted_value = (value - mean_intensity) * T1_contrast + mean_intensity;

            // Normalize T1 voxel value for grayscale color
            float gray = std::clamp((adjusted_value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

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

                // Adjust the voxel value using contrast factor
                float adjusted_value_f = (value_f - mean_intensity_f) * fMRI_contrast + mean_intensity_f;

                // Normalize fMRI value for visualization
                float intensity = std::clamp((adjusted_value_f - minValue_f) / (maxValue_f - minValue_f), 0.0f, 1.0f);

                // Update overlay intensity based on the fMRI value
                overlayIntensity = intensity * fMRI_overlayOpacity;
            }

            // Combine the T1 and fMRI colors
            r = r * (1.0f - overlayIntensity) + overlayIntensity;

            for (int ROI_num = 0; ROI_num < ROI_vector.size(); ROI_num++) {
                if(!ROI_toggle_states[ROI_num]) continue;

                NIBR::Image<bool> &ROI_image = *(ROI_vector[ROI_num]);
                QColor ROI_color = ROI_colors[ROI_num];
                float ROI_opacity = ROI_opacities[ROI_num];
                if (ROI_image.data[idx]) {
                    // Extract RGB components from the ROI color
                    float roi_r = ROI_color.redF();   // Normalize to range [0, 1]
                    float roi_g = ROI_color.greenF(); // Normalize to range [0, 1]
                    float roi_b = ROI_color.blueF();  // Normalize to range [0, 1]
            
                    // Blend the current color with the ROI color using alpha blending (ROI_opacity)
                    r = r * (1.0f - ROI_opacity) + roi_r * ROI_opacity;
                    g = g * (1.0f - ROI_opacity) + roi_g * ROI_opacity;
                    b = b * (1.0f - ROI_opacity) + roi_b * ROI_opacity;
            
                    // Adjust the alpha value (transparency) based on ROI_opacity
                    a = a * (1.0f - ROI_opacity) + ROI_opacity;
                }
            }

            // Set the final color with alpha blending
            glColor4f(r, g, b, a);

            // Adjust coordinates to center around (0,0)
            float x0 = (x * pixDimX_T1) - (imgWidth_coronal / 2.0f);
            float y0 = (z * pixDimZ_T1) - (imgHeight_coronal / 2.0f);
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
    lineX = (xCrosshair * pixDimX_T1) - (imgWidth_coronal / 2.0f) + (pixDimX_T1 / 2.0f);
    
    // Adjust z index for orientation
    yCrosshair = k;
    if (T1_orientation[2] == "I") yCrosshair = dimZ - 1 - k;
    
    // Compute lineY in physical coordinates
    lineY = (yCrosshair * pixDimZ_T1) - (imgHeight_coronal / 2.0f) + (pixDimZ_T1 / 2.0f);
    
    // Vertical line at x = lineX
    glVertex2f(lineX, -imgHeight_coronal / 2.0f);
    glVertex2f(lineX, imgHeight_coronal / 2.0f);
    
    // Horizontal line at z = lineY
    glVertex2f(-imgWidth_coronal / 2.0f, lineY);
    glVertex2f(imgWidth_coronal / 2.0f, lineY);
    
    glEnd();

    // --------------------
    // Draw Sagittal Slice
    // --------------------
    glViewport(sagittalViewport.x, sagittalViewport.y, sagittalViewport.width, widgetHeight);

    // Calculate the physical dimensions of the image in the sagittal plane
    imgWidth_sagittal = dimY * pixDimY_T1;
    imgHeight_sagittal = dimZ * pixDimZ_T1;
    
    // Use the panOffset and zoomFactor
    setupOrtho(imgWidth_sagittal, imgHeight_sagittal, -panOffset.y(), panOffset.z(),
               sagittalViewport.width, sagittalViewport.height,
               leftSagittal, rightSagittal, bottomSagittal, topSagittal);

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

            // Adjust the voxel value using contrast factor
            float adjusted_value = (value - mean_intensity) * T1_contrast + mean_intensity;

            // Normalize T1 voxel value for grayscale color
            float gray = std::clamp((adjusted_value - minValue) / (maxValue - minValue), 0.0f, 1.0f);

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

                // Adjust the voxel value using contrast factor
                float adjusted_value_f = (value_f - mean_intensity_f) * fMRI_contrast + mean_intensity_f;

                // Normalize fMRI value for visualization
                float intensity = std::clamp((adjusted_value_f - minValue_f) / (maxValue_f - minValue_f), 0.0f, 1.0f);

                // Update overlay intensity based on the fMRI value
                overlayIntensity = intensity * fMRI_overlayOpacity;
            }

            // Combine the T1 and fMRI colors
            r = r * (1.0f - overlayIntensity) + overlayIntensity;

            for (int ROI_num = 0; ROI_num < ROI_vector.size(); ROI_num++) {
                if(!ROI_toggle_states[ROI_num]) continue;

                NIBR::Image<bool> &ROI_image = *(ROI_vector[ROI_num]);
                QColor ROI_color = ROI_colors[ROI_num];
                float ROI_opacity = ROI_opacities[ROI_num];
                if (ROI_image.data[idx]) {
                    // Extract RGB components from the ROI color
                    float roi_r = ROI_color.redF();   // Normalize to range [0, 1]
                    float roi_g = ROI_color.greenF(); // Normalize to range [0, 1]
                    float roi_b = ROI_color.blueF();  // Normalize to range [0, 1]
            
                    // Blend the current color with the ROI color using alpha blending (ROI_opacity)
                    r = r * (1.0f - ROI_opacity) + roi_r * ROI_opacity;
                    g = g * (1.0f - ROI_opacity) + roi_g * ROI_opacity;
                    b = b * (1.0f - ROI_opacity) + roi_b * ROI_opacity;
            
                    // Adjust the alpha value (transparency) based on ROI_opacity
                    a = a * (1.0f - ROI_opacity) + ROI_opacity;
                }
            }

            // Set the final color with alpha blending
            glColor4f(r, g, b, a);

            // Adjust coordinates to center around (0,0)
            float x0 = (y * pixDimY_T1) - (imgWidth_sagittal / 2.0f);
            float y0 = (z * pixDimZ_T1) - (imgHeight_sagittal / 2.0f);
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
    lineX = (xCrosshair * pixDimY_T1) - (imgWidth_sagittal / 2.0f) + (pixDimY_T1 / 2.0f);
    
    // Adjust z index for orientation
    yCrosshair = k;
    if (T1_orientation[2] == "I") yCrosshair = dimZ - 1 - k;
    
    // Compute lineY in physical coordinates
    lineY = (yCrosshair * pixDimZ_T1) - (imgHeight_sagittal / 2.0f) + (pixDimZ_T1 / 2.0f);
    
    // Vertical line at x = lineX
    glVertex2f(lineX, -imgHeight_sagittal / 2.0f);
    glVertex2f(lineX, imgHeight_sagittal / 2.0f);
    
    // Horizontal line at z = lineY
    glVertex2f(-imgWidth_sagittal / 2.0f, lineY);
    glVertex2f(imgWidth_sagittal / 2.0f, lineY);
    
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
    painter.setRenderHint(QPainter::Antialiasing);

    // Set text properties
    painter.setPen(Qt::red);
    painter.setFont(QFont("Arial", 10));

    // Draw direction labels for each viewport
    drawDirectionLabels(painter, axialViewport, "axial");
    drawDirectionLabels(painter, coronalViewport, "coronal");
    drawDirectionLabels(painter, sagittalViewport, "sagittal");

    // Set text properties
    painter.setPen(Qt::yellow);

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

void OpenGLMRIWidget::drawDirectionLabels(QPainter &painter, const Viewport &viewport, const QString &sliceType)
{
    // Translate the painter to the viewport origin
    painter.save();
    painter.translate(viewport.x, viewport.y);

    // Get the size of the viewport
    int vpWidth = viewport.width;
    int vpHeight = viewport.height;

    // Determine the labels based on slice type and orientation
    QString labelTop, labelBottom, labelLeft, labelRight;

    if (sliceType == "axial") {
        labelTop = T1_orientation[1] == "A" ? "A" : "P";
        labelBottom = T1_orientation[1] == "A" ? "P" : "A";
        labelLeft = T1_orientation[0] == "R" ? "R" : "L";
        labelRight = T1_orientation[0] == "R" ? "L" : "R";
    } else if (sliceType == "coronal") {
        labelTop = T1_orientation[2] == "S" ? "S" : "I";
        labelBottom = T1_orientation[2] == "S" ? "I" : "S";
        labelLeft = T1_orientation[0] == "R" ? "R" : "L";
        labelRight = T1_orientation[0] == "R" ? "L" : "R";
    } else if (sliceType == "sagittal") {
        labelTop = T1_orientation[2] == "S" ? "S" : "I";
        labelBottom = T1_orientation[2] == "S" ? "I" : "S";
        labelLeft = T1_orientation[1] == "A" ? "A" : "P";
        labelRight = T1_orientation[1] == "A" ? "P" : "A";
    }

    // Set up font metrics
    QFontMetrics fm(painter.font());
    int textWidth = fm.horizontalAdvance("M"); // Assume all labels are single letters
    int textHeight = fm.height();

    int margin = 10; // Margin from the edge

    // Positions for labels
    // Top label
    int topX = (vpWidth - textWidth) / 2;
    int topY = margin + textHeight;

    // Bottom label
    int bottomX = (vpWidth - textWidth) / 2;
    int bottomY = vpHeight - margin;

    // Left label
    int leftX = margin;
    int leftY = (vpHeight + textHeight) / 2;

    // Right label
    int rightX = vpWidth - margin - textWidth;
    int rightY = (vpHeight + textHeight) / 2;

    // Draw labels
    painter.drawText(topX, topY, labelTop);
    painter.drawText(bottomX, bottomY, labelBottom);
    painter.drawText(leftX, leftY, labelLeft);
    painter.drawText(rightX, rightY, labelRight);

    painter.restore();
}


void OpenGLMRIWidget::handleAxialClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    convertAxialScreenToImage(mousePos, viewportWidth, viewportHeight, i, j);
    
    edit_ROIs();

    update(); // Trigger a repaint
}

void OpenGLMRIWidget::handleCoronalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    convertCoronalScreenToImage(mousePos, viewportWidth, viewportHeight, i, k);

    edit_ROIs();

    update(); // Trigger a repaint
}

void OpenGLMRIWidget::handleSagittalClick(const QPoint& mousePos, int viewportWidth, int viewportHeight)
{
    convertSagittalScreenToImage(mousePos, viewportWidth, viewportHeight, j, k);

    edit_ROIs();

    update(); // Trigger a repaint
}

void OpenGLMRIWidget::edit_ROIs() {
    if (editMode) {

        switch (editorMode) {
            case BRUSH_SQUARE:
                paint_square_ROI(brush_size);
                break;

            case BRUSH_CIRCLE:
                paint_circle_ROI(brush_size);
                break;

            case RECTANGLE:
                break;
        }
    }
    redo_stack.clear();
}

void OpenGLMRIWidget::mousePressEvent(QMouseEvent *event)
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
        int widgetWidth = width();
        int viewportWidth = widgetWidth / 3;
        int viewportIndex = event->pos().x() / viewportWidth;

        switch (viewportIndex) {
            case 0:
                currentSliceType = AXIAL;
                break;
            case 1:
                currentSliceType = CORONAL;
                break;
            case 2:
                currentSliceType = SAGITTAL;
                break;
        }

        if (editMode && editorMode == RECTANGLE) {
            // Start drawing rectangle
            rectStartPoint = event->pos();
            rectStartImageCoords = screenToImageCoordinates(rectStartPoint);

            isDrawingRect = true;
        } else {
            // Existing code for handling clicks
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
}


void OpenGLMRIWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (T1_image.data == nullptr || T1_image.voxCnt == 0) {
        // No image loaded; ignore the event
        return;
    }

    if (!mousePressed) {
        return;
    }

    if (event->modifiers() & Qt::ControlModifier) {
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

    } else if ((pressedButton == Qt::RightButton) && (event->modifiers() & Qt::ShiftModifier)) {
        // Handle fMRI contrast adjustment
        QPoint currentPoint = event->pos();
        int deltaY = currentPoint.y() - lastContrastPoint.y(); // Vertical movement
        lastContrastPoint = currentPoint;

        // Adjust the fMRI contrast factor
        float sensitivity = 0.01f; // Adjust this value to control the speed of contrast change

        fMRI_contrast += -deltaY * sensitivity;

        // Clamp the contrast factor to reasonable limits
        fMRI_contrast = std::clamp(fMRI_contrast, 0.1f, 10.0f);

        update(); // Repaint the image with updated contrast

    } else if (pressedButton == Qt::RightButton) {
        // Handle T1 contrast adjustment
        QPoint currentPoint = event->pos();
        int deltaY = currentPoint.y() - lastContrastPoint.y(); // Vertical movement
        lastContrastPoint = currentPoint;

        // Adjust the T1 contrast factor
        float sensitivity = 0.01f; // Adjust this value to control the speed of contrast change

        T1_contrast += -deltaY * sensitivity;

        // Clamp the contrast factor to reasonable limits
        T1_contrast = std::clamp(T1_contrast, 0.1f, 10.0f);

        update(); // Repaint the image with updated contrast
    }
}

void OpenGLMRIWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if (T1_image.data == nullptr || T1_image.voxCnt == 0) {
        // No image loaded; ignore the event
        return;
    }

    Q_UNUSED(event);
    
    if (pressedButton == Qt::LeftButton && editMode && editorMode == RECTANGLE && isDrawingRect) {
        // Finish drawing rectangle
        rectCurrentPoint = event->pos();

        // Convert the mouse position to image coordinates
        QPoint imageCoords = screenToImageCoordinates(rectCurrentPoint);
        rectEndImageCoords = imageCoords;

        isDrawingRect = false;

        // Update the ROI between the two points
        applyRectangleToROI(rectStartImageCoords, rectEndImageCoords);

        update(); // Trigger a repaint
    }

    if (editMode && undo_stack_temp.size() > 0) {
        std::vector<int64_t> changes(undo_stack_temp.begin(), undo_stack_temp.end());
        undo_stack.push_back(changes);
        undo_stack_temp.clear();
    }

    mousePressed = false;
}

void OpenGLMRIWidget::handleAxialScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120; // Each notch is 15 degrees; 120 units per notch
    k = std::clamp(k + delta, 0, static_cast<int>(T1_image.imgDims[2] - 1)); // Adjust axial index (z-axis)
    update();
}

void OpenGLMRIWidget::handleCoronalScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120;
    j = std::clamp(j + delta, 0, static_cast<int>(T1_image.imgDims[1] - 1)); // Adjust coronal index (y-axis)
    update();
}

void OpenGLMRIWidget::handleSagittalScroll(QWheelEvent *event)
{
    int delta = event->angleDelta().y() / 120;
    i = std::clamp(i + delta, 0, static_cast<int>(T1_image.imgDims[0] - 1)); // Adjust sagittal index (x-axis)
    update();
}

void OpenGLMRIWidget::wheelEvent(QWheelEvent *event)
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

void OpenGLMRIWidget::keyPressEvent(QKeyEvent *event)
{
    if ((event->key() == Qt::Key_Z) && (event->modifiers() & Qt::ControlModifier)) {
        // Handle Ctrl+Z (Undo)
        undoButton_clicked(); // Call your undo function
    }
    else if ((event->key() == Qt::Key_Y) && (event->modifiers() & Qt::ControlModifier)) {
        // Handle Ctrl+Y (Redo)
        redoButton_clicked(); // Call your redo function
    }
    else if (event->key() == Qt::Key_R) {
        // Reset zoom and pan
        zoomFactor = 1.0f;
        panOffset = QVector3D(0.0f, 0.0f, 0.0f);
        T1_contrast = 1.0f;
        fMRI_contrast = 1.0f;
        update();
    }
    else {
        // Call the base class implementation for unhandled keys
        QWidget::keyPressEvent(event);
    }
}

void OpenGLMRIWidget::applyRectangleToROI(const QPoint& startImageCoords, const QPoint& endImageCoords)
{
    int imgXStart = startImageCoords.x();
    int imgYStart = startImageCoords.y();
    int imgXEnd = endImageCoords.x();
    int imgYEnd = endImageCoords.y();

    // Ensure start indices are less than end indices
    int imgXMin = std::min(imgXStart, imgXEnd);
    int imgXMax = std::max(imgXStart, imgXEnd);
    int imgYMin = std::min(imgYStart, imgYEnd);
    int imgYMax = std::max(imgYStart, imgYEnd);

    // Get the current slice index (k) for axial, (j) for coronal, (i) for sagittal
    int fixedIndex;
    switch (currentSliceType) {
        case AXIAL:
            fixedIndex = k; // Current axial slice index
            break;
        case CORONAL:
            fixedIndex = j; // Current coronal slice index
            break;
        case SAGITTAL:
            fixedIndex = i; // Current sagittal slice index
            break;
    }

    // Update the ROI
    NIBR::Image<bool> &ROI_image = *(ROI_vector[target_ROI]);

    for (int x = imgXMin; x <= imgXMax; x++) {
        for (int y = imgYMin; y <= imgYMax; y++) {
            int idxi, idxj, idxk;

            switch (currentSliceType) {
                case AXIAL:
                    idxi = x;
                    idxj = y;
                    idxk = fixedIndex;
                    break;
                case CORONAL:
                    idxi = x;
                    idxj = fixedIndex;
                    idxk = y;
                    break;
                case SAGITTAL:
                    idxi = fixedIndex;
                    idxj = x;
                    idxk = y;
                    break;
            }

            int64_t target_index = ROI_image.sub2ind(idxi, idxj, idxk);
            if (ROI_image.data[target_index] == isMarking) continue;
            else {
                ROI_image.data[target_index] = isMarking;
                undo_stack_temp.insert(target_index);
            }
        }
    }
    redo_stack.clear();
}

QPoint OpenGLMRIWidget::screenToImageCoordinates(const QPoint& screenPoint)
{
    // Determine which viewport the point is in
    int widgetWidth = width();
    int viewportWidth = widgetWidth / 3;
    int viewportIndex = screenPoint.x() / viewportWidth;

    int imgX, imgY;
    switch (viewportIndex) {
        case 0:
            // Axial slice
            convertAxialScreenToImage(screenPoint, viewportWidth, height(), imgX, imgY);
            return QPoint(imgX, imgY);
        case 1:
            // Coronal slice
            convertCoronalScreenToImage(screenPoint, viewportWidth, height(), imgX, imgY);
            return QPoint(imgX, imgY);
        case 2:
            // Sagittal slice
            convertSagittalScreenToImage(screenPoint, viewportWidth, height(), imgX, imgY);
            return QPoint(imgX, imgY);
        default:
            return QPoint(0, 0); // Default case
    }
}

void OpenGLMRIWidget::convertAxialScreenToImage(const QPoint& screenPoint, int viewportWidth, int viewportHeight, int& imgi, int& imgj)
{
    // Dimensions of the Axial slice
    int dimX = T1_image.imgDims[0];
    int dimY = T1_image.imgDims[1];

    // Calculate the position within the axial viewport
    int xInViewport = screenPoint.x(); // Since axialViewport.x = 0
    int yInViewport = screenPoint.y();

    // Flip y-coordinate because Qt's y=0 is at the top, OpenGL's y=0 is at the bottom
    int yInViewportFlipped = viewportHeight - yInViewport;

    // Convert viewport coordinates to Normalized Device Coordinates (NDC)
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewportFlipped) / viewportHeight) * 2.0f - 1.0f;

    // Map NDC to world coordinates using stored projection parameters
    float worldX = ((ndcX + 1.0f) / 2.0f) * (rightAxial - leftAxial) + leftAxial;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (topAxial - bottomAxial) + bottomAxial;

    // Convert world coordinates to image coordinates
    float imgXf = (worldX + imgWidth_axial / 2.0f) / T1_pixDims[0];
    float imgYf = (worldY + imgHeight_axial / 2.0f) / T1_pixDims[1];

    int imgX = static_cast<int>(imgXf);
    int imgY = static_cast<int>(imgYf);

    // Adjust for image orientation
    if (T1_orientation[0] == "R") imgX = dimX - 1 - imgX;
    if (T1_orientation[1] == "P") imgY = dimY - 1 - imgY;

    // Clamp the indices to valid ranges
    imgi = std::clamp(imgX, 0, dimX - 1);
    imgj = std::clamp(imgY, 0, dimY - 1);
}

void OpenGLMRIWidget::convertCoronalScreenToImage(const QPoint& screenPoint, int viewportWidth, int viewportHeight, int& imgi, int& imgk)
{
    // Dimensions of the Coronal slice
    int dimX = T1_image.imgDims[0];
    int dimZ = T1_image.imgDims[2];

    // Calculate the position within the coronal viewport
    int xInViewport = screenPoint.x() - viewportWidth; // Since coronal viewport starts at viewportWidth
    int yInViewport = screenPoint.y();

    // Flip y-coordinate
    int yInViewportFlipped = viewportHeight - yInViewport;

    // Convert viewport coordinates to NDC
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewportFlipped) / viewportHeight) * 2.0f - 1.0f;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (rightCoronal - leftCoronal) + leftCoronal;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (topCoronal - bottomCoronal) + bottomCoronal;

    // Convert world coordinates to image coordinates
    float imgXf = (worldX + imgWidth_coronal / 2.0f) / T1_pixDims[0];
    float imgZf = (worldY + imgHeight_coronal / 2.0f) / T1_pixDims[2];

    int imgX = static_cast<int>(imgXf);
    int imgZ = static_cast<int>(imgZf);

    // Adjust for image orientation
    if (T1_orientation[0] == "R") imgX = dimX - 1 - imgX;
    if (T1_orientation[2] == "I") imgZ = dimZ - 1 - imgZ;

    // Clamp the indices
    imgi = std::clamp(imgX, 0, dimX - 1);
    imgk = std::clamp(imgZ, 0, dimZ - 1);
}

void OpenGLMRIWidget::convertSagittalScreenToImage(const QPoint& screenPoint, int viewportWidth, int viewportHeight, int& imgj, int& imgk)
{
    // Dimensions of the Coronal slice
    int dimY = T1_image.imgDims[1];
    int dimZ = T1_image.imgDims[2];

    // Calculate the position within the coronal viewport
    int xInViewport = screenPoint.x() - 2 * viewportWidth;
    int yInViewport = screenPoint.y();

    // Flip y-coordinate
    int yInViewportFlipped = viewportHeight - yInViewport;

    // Convert viewport coordinates to NDC
    float ndcX = (static_cast<float>(xInViewport) / viewportWidth) * 2.0f - 1.0f;
    float ndcY = (static_cast<float>(yInViewportFlipped) / viewportHeight) * 2.0f - 1.0f;

    // Map NDC to world coordinates
    float worldX = ((ndcX + 1.0f) / 2.0f) * (rightSagittal - leftSagittal) + leftSagittal;
    float worldY = ((ndcY + 1.0f) / 2.0f) * (topSagittal - bottomSagittal) + bottomSagittal;

    // Convert world coordinates to image coordinates
    float imgYf = (worldX + imgWidth_sagittal / 2.0f) / T1_pixDims[1];
    float imgZf = (worldY + imgHeight_sagittal / 2.0f) / T1_pixDims[2];

    int imgY = static_cast<int>(imgYf);
    int imgZ = static_cast<int>(imgZf);

    // Adjust for image orientation
    if (T1_orientation[1] == "A") imgY = dimY - 1 - imgY;
    if (T1_orientation[2] == "I") imgZ = dimZ - 1 - imgZ;

    // Clamp the indices
    imgj = std::clamp(imgY, 0, dimY - 1);
    imgk = std::clamp(imgZ, 0, dimZ - 1);
}

void OpenGLMRIWidget::updateGraph()
{
    // This method should ideally handle fetching new data and triggering a redraw
    update();  // Request a re-draw
}













// ----------------------------
// ROI functions
// ----------------------------
void OpenGLMRIWidget::addButton_clicked() 
{
    // Create a new shared_ptr for the image
    auto newImage = std::make_shared<NIBR::Image<bool>>();
    newImage->createFromTemplate(T1_image, true);

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
    ROI_vector.push_back(newImage); // Add the shared_ptr to the vector
    ROI_toggle_states.push_back(false);
    ROI_colors.push_back(QColor(0, 255, 0));
    ROI_opacities.push_back(0.5f);

    emit ROI_update(ROI_names, ROI_toggle_states);
    emit ROIListChanged(ROI_names);  // Add this
    reset_undo_stacks();
    update();

    ROI_fMRI_means.conservativeResize(ROI_fMRI_means.size() + 1);  // Add new element
    ROI_fMRI_means(ROI_fMRI_means.size() - 1) = 0.0;  // Set new element to 0
    if (fMRI_image_resampled.voxCnt > 0) {
        calculateROIMeans();  // Recalculate means if fMRI data exists
    }
}



void OpenGLMRIWidget::loadButton_clicked(const QString& filePath) 
{
    std::string filePath_std = filePath.toStdString();

    // Create a new shared_ptr for the image
    auto newImage = std::make_shared<NIBR::Image<bool>>(filePath_std);
    newImage->read();

    ROI_vector.push_back(newImage); // Add the shared_ptr to the vector

    // Extract the file name from the file path
    size_t lastSlash = filePath_std.find_last_of("/\\");
    size_t lastDot = filePath_std.find_last_of(".");
    std::string fileName = filePath_std.substr(lastSlash + 1, lastDot - lastSlash - 1);

    ROI_names.push_back(fileName);
    ROI_toggle_states.push_back(false);
    ROI_colors.push_back(QColor(0, 255, 0));
    ROI_opacities.push_back(0.5f);

    emit ROI_update(ROI_names, ROI_toggle_states);
    emit ROIListChanged(ROI_names);  // Add this
    reset_undo_stacks();
    update();
}


void OpenGLMRIWidget::saveButton_clicked() 
{
    if (target_ROI >= 0 && target_ROI < ROI_vector.size()) {
        ROI_vector[target_ROI]->write(ROI_names[target_ROI] + ".nii.gz");
    } else {
        // Handle invalid target_ROI
        QMessageBox::warning(this, "Warning", "No ROI selected to save.");
    }
}


void OpenGLMRIWidget::deleteButton_clicked() 
{
    if (target_ROI >= 0 && target_ROI < ROI_vector.size()) {
        ROI_names.erase(ROI_names.begin() + target_ROI);
        ROI_vector.erase(ROI_vector.begin() + target_ROI);
        ROI_toggle_states.erase(ROI_toggle_states.begin() + target_ROI);
        ROI_colors.erase(ROI_colors.begin() + target_ROI);
        ROI_opacities.erase(ROI_opacities.begin() + target_ROI);

        // Adjust target_ROI if necessary
        if (target_ROI >= ROI_vector.size()) {
            target_ROI = ROI_vector.empty() ? -1 : ROI_vector.size() - 1;
        }

        if (ROI_vector.empty()) {
            editMode = false;
        }

        // Debugging output
        for (size_t i = 0; i < ROI_vector.size(); ++i) {
            print_debug("ROI " + std::to_string(i) + ": " + ROI_names[i] + " toggle state: " + std::to_string(ROI_toggle_states[i]));
        }

        emit ROI_update(ROI_names, ROI_toggle_states);
        emit ROIListChanged(ROI_names);  // Add this
        reset_undo_stacks();
        update();

        if (target_ROI >= 0 && target_ROI < ROI_fMRI_means.size()) {
            // Create new vector without the element at target_ROI
            Eigen::VectorXd newMeans(ROI_fMRI_means.size() - 1);
            if (target_ROI > 0) {
                newMeans.segment(0, target_ROI) = ROI_fMRI_means.segment(0, target_ROI);
            }
            if (target_ROI < ROI_fMRI_means.size() - 1) {
                newMeans.segment(target_ROI, ROI_fMRI_means.size() - target_ROI - 1) = 
                    ROI_fMRI_means.segment(target_ROI + 1, ROI_fMRI_means.size() - target_ROI - 1);
            }
            ROI_fMRI_means = newMeans;
        }
    } else {
        // Handle invalid target_ROI
        QMessageBox::warning(this, "Warning", "No ROI selected to delete.");
    }
}


void OpenGLMRIWidget::onSliderValueChanged(int value) 
{
    ROI_opacities[target_ROI] = value / 100.0f;
    update();
}

void OpenGLMRIWidget::colorButton_clicked(QColor color) 
{
    ROI_colors[target_ROI] = color;
    update();
}

void OpenGLMRIWidget::editButton_toggled(bool checked)
{
    if (ROI_vector.size() > 0) {
        editMode = checked;
    }
}

void OpenGLMRIWidget::undoButton_clicked()
{
    if (!undo_stack.empty()) {
        // Get the last change from the undo stack
        std::vector<int64_t> last_change = undo_stack.back();
        undo_stack.pop_back(); // Remove it from the undo stack

        // Revert the changes
        for (int64_t index : last_change) {
            revert_ROI_at_index(index); // Revert the change at index
        }

        // Push the change onto the redo stack
        redo_stack.push_back(last_change);
    }

    update();
}

void OpenGLMRIWidget::redoButton_clicked()
{
    if (!redo_stack.empty()) {
        // Get the last undone change from the redo stack
        std::vector<int64_t> last_redo = redo_stack.back();
        redo_stack.pop_back(); // Remove it from the redo stack

        // Reapply the changes
        for (int64_t index : last_redo) {
            revert_ROI_at_index(index); // Reapply the change at index
        }

        // Push the change back onto the undo stack
        undo_stack.push_back(last_redo);
    }

    update();
}

void OpenGLMRIWidget::revert_ROI_at_index(int64_t image_index) {
    NIBR::Image<bool> &ROI_image = *(ROI_vector[target_ROI]);
    if (ROI_image.data[image_index] == 1) ROI_image.data[image_index] = 0;
    else ROI_image.data[image_index] = 1;
}

void OpenGLMRIWidget::paint_square_ROI(int brush_size) {
    // Get image dimensions for bounds checking
    int dimX = T1_image.imgDims[0];
    int dimY = T1_image.imgDims[1];
    int dimZ = T1_image.imgDims[2];

    switch (currentSliceType) {
        case AXIAL:
            for (int i_temp = i - brush_size; i_temp <= i + brush_size; i_temp++) {
                for (int j_temp = j - brush_size; j_temp <= j + brush_size; j_temp++) {
                    // Check bounds to prevent accessing invalid indices
                    if (i_temp >= 0 && i_temp < dimX &&
                        j_temp >= 0 && j_temp < dimY) {

                        int64_t target_index = T1_image.sub2ind(i_temp, j_temp, k);
                        if (!(T1_image.data[target_index] == isMarking)) {
                            set_ROI_at_index(target_index, isMarking);
                            undo_stack_temp.insert(target_index);
                        }
                    }
                }
            }
            break;

        case CORONAL:
            for (int i_temp = i - brush_size; i_temp <= i + brush_size; i_temp++) {
                for (int k_temp = k - brush_size; k_temp <= k + brush_size; k_temp++) {
                    // Check bounds
                    if (i_temp >= 0 && i_temp < dimX &&
                        k_temp >= 0 && k_temp < dimZ) {

                        int64_t target_index = T1_image.sub2ind(i_temp, j, k_temp);
                        if (!(T1_image.data[target_index] == isMarking)) {
                            set_ROI_at_index(target_index, isMarking);
                            undo_stack_temp.insert(target_index);
                        }
                    }
                }
            }
            break;

        case SAGITTAL:
            for (int j_temp = j - brush_size; j_temp <= j + brush_size; j_temp++) {
                for (int k_temp = k - brush_size; k_temp <= k + brush_size; k_temp++) {
                    // Check bounds
                    if (j_temp >= 0 && j_temp < dimY &&
                        k_temp >= 0 && k_temp < dimZ) {

                        int64_t target_index = T1_image.sub2ind(i, j_temp, k_temp);
                        if (!(T1_image.data[target_index] == isMarking)) {
                            set_ROI_at_index(target_index, isMarking);
                            undo_stack_temp.insert(target_index);
                        }
                    }
                }
            }
            break;
    }
}

void OpenGLMRIWidget::paint_circle_ROI(int brush_size) {
    switch (currentSliceType) {
        case AXIAL:
            for (int i_temp = i - brush_size; i_temp <= i + brush_size; i_temp++) {
                for (int j_temp = j - brush_size; j_temp <= j + brush_size; j_temp++) {
                    // Calculate the squared distance from the center
                    int di = i_temp - i;
                    int dj = j_temp - j;
                    if (di*di + dj*dj <= brush_size*brush_size) {
                        // Check bounds to prevent accessing invalid indices
                        if (i_temp >= 0 && i_temp < T1_image.imgDims[0] &&
                            j_temp >= 0 && j_temp < T1_image.imgDims[1]) {

                            int64_t target_index = T1_image.sub2ind(i_temp, j_temp, k);
                            if (!(T1_image.data[target_index] == isMarking)) {
                                set_ROI_at_index(target_index, isMarking);
                                undo_stack_temp.insert(target_index);
                            }
                        }
                    }
                }
            }
            break;

        case CORONAL:
            for (int i_temp = i - brush_size; i_temp <= i + brush_size; i_temp++) {
                for (int k_temp = k - brush_size; k_temp <= k + brush_size; k_temp++) {
                    // Calculate the squared distance from the center
                    int di = i_temp - i;
                    int dk = k_temp - k;
                    if (di*di + dk*dk <= brush_size*brush_size) {
                        // Check bounds
                        if (i_temp >= 0 && i_temp < T1_image.imgDims[0] &&
                            k_temp >= 0 && k_temp < T1_image.imgDims[2]) {

                            int64_t target_index = T1_image.sub2ind(i_temp, j, k_temp);
                            if (!(T1_image.data[target_index] == isMarking)) {
                                set_ROI_at_index(target_index, isMarking);
                                undo_stack_temp.insert(target_index);
                            }
                        }
                    }
                }
            }
            break;

        case SAGITTAL:
            for (int j_temp = j - brush_size; j_temp <= j + brush_size; j_temp++) {
                for (int k_temp = k - brush_size; k_temp <= k + brush_size; k_temp++) {
                    // Calculate the squared distance from the center
                    int dj = j_temp - j;
                    int dk = k_temp - k;
                    if (dj*dj + dk*dk <= brush_size*brush_size) {
                        // Check bounds
                        if (j_temp >= 0 && j_temp < T1_image.imgDims[1] &&
                            k_temp >= 0 && k_temp < T1_image.imgDims[2]) {

                            int64_t target_index = T1_image.sub2ind(i, j_temp, k_temp);
                            if (!(T1_image.data[target_index] == isMarking)) {
                                set_ROI_at_index(target_index, isMarking);
                                undo_stack_temp.insert(target_index);
                            }
                        }
                    }
                }
            }
            break;
    }
}

void OpenGLMRIWidget::set_ROI_at_index(int64_t image_index, bool ROI_value) {
        NIBR::Image<bool> &ROI_image = *(ROI_vector[target_ROI]);
        ROI_image.data[image_index] = ROI_value;
}

void OpenGLMRIWidget::aboveButton_clicked() {
    if (!undo_stack.empty() && !ROI_vector.empty()) {
        // Get the last modification
        std::vector<int64_t> last_change = undo_stack.back();

        // Get image dimensions
        NIBR::Image<bool>& ROI_image = *(ROI_vector[target_ROI]); // Assuming all images have the same dimensions
        int dimX = ROI_image.imgDims[0];
        int dimY = ROI_image.imgDims[1];
        int dimZ = ROI_image.imgDims[2];

        // Determine the incremented index and maximum allowed index
        int incrementedIndex;
        int maxIndex;

        // Calculate the incremented index based on the current slice type
        switch (currentSliceType) {
            case AXIAL:
                k = k + 1;
                incrementedIndex = k;
                maxIndex = dimZ - 1;
                break;
            case CORONAL:
                j = j + 1;
                incrementedIndex = j;
                maxIndex = dimY - 1;
                break;
            case SAGITTAL:
                i = i + 1;
                incrementedIndex = i;
                maxIndex = dimX - 1;
                break;
            default:
                std::cout << "Invalid slice type: " << currentSliceType << std::endl;
                return;
        }

        // Check if we can move "above"
        if (incrementedIndex > maxIndex) {
            std::cout << "Already at the top slice. Cannot copy to above slice." << std::endl;
            return;
        }

        // Set to store unique indices for the above slice
        std::unordered_set<int64_t> above_indices_set;

        // For each modified index, compute the corresponding index in the above slice
        for (int64_t index : last_change) {
            int64_t idxi, idxj, idxk;

            // Convert linear index to (idxi, idxj, idxk)
            ROI_image.ind2sub(index, idxi, idxj, idxk);

            // Adjust the appropriate index
            switch (currentSliceType) {
                case AXIAL:
                    idxk = incrementedIndex;
                    break;
                case CORONAL:
                    idxj = incrementedIndex;
                    break;
                case SAGITTAL:
                    idxi = incrementedIndex;
                    break;
            }

            // Check bounds
            if (idxi < 0 || idxi >= dimX || idxj < 0 || idxj >= dimY || idxk < 0 || idxk >= dimZ)
                continue;

            // Compute the new linear index
            int64_t above_index = ROI_image.sub2ind(idxi, idxj, idxk);

            // Insert into the set to ensure uniqueness
            above_indices_set.insert(above_index);

            // Apply the modification
            set_ROI_at_index(above_index, true);
        }

        // Convert the set to a vector and push onto the undo stack
        std::vector<int64_t> above_indices(above_indices_set.begin(), above_indices_set.end());
        undo_stack.push_back(above_indices);

        // Clear the redo stack since we have a new modification
        redo_stack.clear();

        // Optionally, update the current slice index if needed
        // currentSlice = incrementedIndex;

        // Refresh the display if necessary
        update();
    } else {
        std::cout << "No modifications to copy from." << std::endl;
    }
}

void OpenGLMRIWidget::belowButton_clicked() {
    if (!undo_stack.empty() && !ROI_vector.empty()) {
        // Get the last modification
        std::vector<int64_t> last_change = undo_stack.back();

        // Get image dimensions
        NIBR::Image<bool>& ROI_image = *(ROI_vector[target_ROI]); // Assuming all images have the same dimensions
        int dimX = ROI_image.imgDims[0];
        int dimY = ROI_image.imgDims[1];
        int dimZ = ROI_image.imgDims[2];

        // Determine the incremented index and maximum allowed index
        int incrementedIndex;

        // Calculate the incremented index based on the current slice type
        switch (currentSliceType) {
            case AXIAL:
                k = k - 1;
                incrementedIndex = k;
                break;
            case CORONAL:
                j = j - 1;
                incrementedIndex = j;
                break;
            case SAGITTAL:
                i = i - 1;
                incrementedIndex = i;
                break;
            default:
                std::cout << "Invalid slice type: " << currentSliceType << std::endl;
                return;
        }

        // Check if we can move "above"
        if (incrementedIndex < 0 ) {
            std::cout << "Already at the bottom slice. Cannot copy to below slice." << std::endl;
            return;
        }

        // Set to store unique indices for the above slice
        std::unordered_set<int64_t> above_indices_set;

        // For each modified index, compute the corresponding index in the above slice
        for (int64_t index : last_change) {
            int64_t idxi, idxj, idxk;

            // Convert linear index to (idxi, idxj, idxk)
            ROI_image.ind2sub(index, idxi, idxj, idxk);

            // Adjust the appropriate index
            switch (currentSliceType) {
                case AXIAL:
                    idxk = incrementedIndex;
                    break;
                case CORONAL:
                    idxj = incrementedIndex;
                    break;
                case SAGITTAL:
                    idxi = incrementedIndex;
                    break;
            }

            // Check bounds
            if (idxi < 0 || idxi >= dimX || idxj < 0 || idxj >= dimY || idxk < 0 || idxk >= dimZ)
                continue;

            // Compute the new linear index
            int64_t above_index = ROI_image.sub2ind(idxi, idxj, idxk);

            // Insert into the set to ensure uniqueness
            above_indices_set.insert(above_index);

            // Apply the modification
            set_ROI_at_index(above_index, true);
        }

        // Convert the set to a vector and push onto the undo stack
        std::vector<int64_t> above_indices(above_indices_set.begin(), above_indices_set.end());
        undo_stack.push_back(above_indices);

        // Clear the redo stack since we have a new modification
        redo_stack.clear();

        // Refresh the display if necessary
        update();
    } else {
        std::cout << "No modifications to copy from." << std::endl;
    }
}

void OpenGLMRIWidget::calculateROIMeans() 
{
    if (fMRI_image_resampled.voxCnt == 0 || ROI_vector.empty() || !data_handler) {
        return;  // No fMRI data, ROIs to process, or data handler
    }

    // Resize the Eigen vector to match the number of ROIs
    ROI_fMRI_means.resize(ROI_vector.size());
    ROI_fMRI_means.setZero();  // Initialize all means to 0

    // For each ROI
    for (size_t roi_idx = 0; roi_idx < ROI_vector.size(); ++roi_idx) {
        NIBR::Image<bool>& ROI_image = *(ROI_vector[roi_idx]);
        double sum = 0.0;  // Changed to double for Eigen compatibility
        int count = 0;

        // For each voxel in the ROI
        for (int64_t i = 0; i < ROI_image.voxCnt; ++i) {
            if (ROI_image.data[i]) {  // If voxel is part of ROI
                sum += fMRI_image_resampled.data[i];
                count++;
            }
        }

        // Calculate mean
        ROI_fMRI_means(roi_idx) = count > 0 ? sum / count : 0.0;
        
        print_debug("ROI " + ROI_names[roi_idx] + " mean fMRI value: " + std::to_string(ROI_fMRI_means(roi_idx)) + " (from " + std::to_string(count) + " voxels)");
    }

    // Set the ROI means in the data handler
    data_handler->setROIMeans(ROI_fMRI_means);
}