#include "mriwindow.h"

mriWindow::mriWindow(QWidget *parent)
    : QMainWindow(parent)
{
    // Create the central widget
    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);

    // Create the main layout
    mainLayout = new QVBoxLayout(centralWidget);

    // Initialize your custom GL widget
    MRIWidget = new OpenGLMRIWidget(centralWidget);
    mainLayout->addWidget(MRIWidget);

    // ----------------------------
    // Create the menu bar
    // ----------------------------
    QMenuBar *menuBar = this->menuBar();

    // Create the Edit menu and add actions
    QMenu *mriMenu = menuBar->addMenu("MRI");
    QAction *T1LoadAction = mriMenu->addAction("Load T1 image");
    QAction *fmriLoadAction = mriMenu->addAction("Load fMRI image");
    QAction *mriSettingsAction = mriMenu->addAction("Tools");

    // Connect actions to slots
    connect(T1LoadAction, &QAction::triggered, this, &mriWindow::MRI_T1_load_clicked);
    connect(fmriLoadAction, &QAction::triggered, this, &mriWindow::fMRI_load_clicked);

    // Set the layout to the central widget
    centralWidget->setLayout(mainLayout);    // Create the dock widget
    QDockWidget *sidePanel = new QDockWidget("ROI controls", this);
    sidePanel->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    // Create the main control widget for the dock
    QWidget *controlWidget = new QWidget();
    QVBoxLayout *mainLayout = new QVBoxLayout(controlWidget);

    // ----------------------------
    // Create the top row of buttons
    // ----------------------------
    QHBoxLayout *buttonLayout = new QHBoxLayout();

    QPushButton *addButton = new QPushButton(QString("add"), controlWidget);
    QPushButton *loadButton = new QPushButton(QString("load"), controlWidget);
    QPushButton *saveButton = new QPushButton(QString("save"), controlWidget);
    QPushButton *deleteButton = new QPushButton(QString("delete"), controlWidget);
    QPushButton *visibleButton = new QPushButton(QString("visible"), controlWidget);

    // Connect buttons to slots
    connect(addButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::addButton_clicked);
    connect(loadButton, &QPushButton::clicked, this, &mriWindow::loadButton_clicked);
    connect(saveButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::saveButton_clicked);
    connect(deleteButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::deleteButton_clicked);
    connect(visibleButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::visibleButton_clicked);

    connect(MRIWidget, &OpenGLMRIWidget::ROI_update, this, &mriWindow::ROI_update);

    addButton->setMinimumWidth(50);
    loadButton->setMinimumWidth(50);
    saveButton->setMinimumWidth(50);
    deleteButton->setMinimumWidth(50);
    visibleButton->setMinimumWidth(50);

    buttonLayout->addWidget(addButton);
    buttonLayout->addWidget(loadButton);
    buttonLayout->addWidget(saveButton);
    buttonLayout->addWidget(deleteButton);
    buttonLayout->addWidget(visibleButton);
    
    mainLayout->addLayout(buttonLayout);

    // ----------------------------
    // Create the list of toggleable items
    // ----------------------------
    toggleList = new QListWidget(controlWidget);
    mainLayout->addWidget(toggleList);

    connect(toggleList, &QListWidget::itemChanged, this, &mriWindow::updateToggleStatus);

    // Create a color picker button
    colorPickerButton = new QPushButton("Color", controlWidget);

    QSlider *slider = new QSlider(Qt::Horizontal, controlWidget);
    slider->setRange(0, 100);
    slider->setValue(50);

    // Create a horizontal layout to place the color picker button and slider side by side
    QHBoxLayout *colorAndSliderLayout = new QHBoxLayout();
    colorAndSliderLayout->addWidget(colorPickerButton); // Add the color picker button
    colorAndSliderLayout->addWidget(slider);            // Add the slider

    // Add the horizontal layout to the main layout
    mainLayout->addLayout(colorAndSliderLayout);

    connect(slider, &QSlider::valueChanged, MRIWidget, &OpenGLMRIWidget::onSliderValueChanged);
    connect(colorPickerButton, &QPushButton::clicked, this, &mriWindow::colorButton_clicked);
    
    // ----------------------------
    // Add Edit, Undo, and Redo buttons
    // ----------------------------
    QHBoxLayout *actionButtonsLayout = new QHBoxLayout();
    QPushButton *editButton = new QPushButton("Edit", controlWidget);
    QPushButton *undoButton = new QPushButton("Undo", controlWidget);
    QPushButton *redoButton = new QPushButton("Redo", controlWidget);
    
    // Make the edit button toggleable
    editButton->setCheckable(true);

    connect(editButton, &QPushButton::toggled, MRIWidget, &OpenGLMRIWidget::editButton_toggled);
    connect(undoButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::undoButton_clicked);
    connect(redoButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::redoButton_clicked);

    connect(editButton, &QPushButton::toggled, [editButton](bool checked) {
        if (checked) { editButton->setStyleSheet("background-color: darkgreen;"); } 
        else         { editButton->setStyleSheet(""); }
    });

    actionButtonsLayout->addWidget(editButton);
    actionButtonsLayout->addWidget(undoButton);
    actionButtonsLayout->addWidget(redoButton);
    mainLayout->addLayout(actionButtonsLayout);

    // ----------------------------
    // Add a frame for editing tools
    // ----------------------------
    QFrame *editFrame = new QFrame(controlWidget);
    editFrame->setFrameShape(QFrame::StyledPanel);
    QVBoxLayout *editFrameLayout = new QVBoxLayout(editFrame);

    // Add Brush and Rectangle selection buttons
    QHBoxLayout *toolButtonsLayout = new QHBoxLayout();
    QPushButton *brushButton = new QPushButton("Brush", editFrame);
    QPushButton *eraseButton = new QPushButton("Erase", editFrame);
    QPushButton *rectangleButton = new QPushButton("Rectangle", editFrame);

    connect(brushButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::brushButton_clicked);
    connect(eraseButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::eraseButton_clicked);
    connect(rectangleButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::rectangleButton_clicked);
    
    toolButtonsLayout->addWidget(brushButton);
    toolButtonsLayout->addWidget(eraseButton);
    toolButtonsLayout->addWidget(rectangleButton);
    editFrameLayout->addLayout(toolButtonsLayout);

    // Add label and Above/Below buttons for "Copy from slice"
    QHBoxLayout *copyLayout = new QHBoxLayout();
    QLabel *copyLabel = new QLabel("Copy from slice:", editFrame);
    QPushButton *aboveButton = new QPushButton("Above", editFrame);
    QPushButton *belowButton = new QPushButton("Below", editFrame);

    connect(aboveButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::aboveButton_clicked);
    connect(belowButton, &QPushButton::clicked, MRIWidget, &OpenGLMRIWidget::belowButton_clicked);

    copyLayout->addWidget(copyLabel);
    copyLayout->addWidget(aboveButton);
    copyLayout->addWidget(belowButton);
    editFrameLayout->addLayout(copyLayout);

    mainLayout->addWidget(editFrame); // Add the frame to the main layout

    // Set the layout for the control widget and assign it to the dock
    controlWidget->setLayout(mainLayout);
    sidePanel->setWidget(controlWidget);

    // Add the dock widget to the main window, docked to the right by default
    addDockWidget(Qt::RightDockWidgetArea, sidePanel);
    sidePanel->setVisible(false);

    // Connect action to toggle the side panel visibility
    connect(mriSettingsAction, &QAction::triggered, [=]() {
        sidePanel->setVisible(!sidePanel->isVisible());
    });
}

mriWindow::~mriWindow()
{
    // No need to delete child widgets, Qt handles it automatically
}

void mriWindow::MRI_T1_load_clicked()
{
    QString filePath = QFileDialog::getOpenFileName(this, "Open MRI Image", "", "NIFTI Images (*.nii *.nii.gz)");
    if (!filePath.isEmpty()) {
        MRIWidget->loadImage_T1(filePath);
    }
}

void mriWindow::fMRI_load_clicked()
{
    QString filePath = QFileDialog::getOpenFileName(this, "Open fMRI Image", "", "NIFTI Images (*.nii *.nii.gz)");
    if (!filePath.isEmpty()) {
        MRIWidget->loadImage_fMRI(filePath);
    }
}

void mriWindow::updateToggleStatus() {
    std::vector<bool> states;
    for (int i = 0; i < toggleList->count(); ++i) {
        QListWidgetItem *item = toggleList->item(i);
        QString name = item->text();
        bool isChecked = (item->checkState() == Qt::Checked);
        states.push_back(isChecked);
    }
    MRIWidget->updateToggleStates(states);
}

void mriWindow::loadButton_clicked() {
    QString filePath = QFileDialog::getOpenFileName(this, "Open MRI Image", "", "NIFTI Images (*.nii *.nii.gz)");
    if (!filePath.isEmpty()) {
        MRIWidget->loadButton_clicked(filePath);
    }
}
void mriWindow::colorButton_clicked() {

    QColor color = QColorDialog::getColor(Qt::white, this, "Select Color");
    if (color.isValid()) {
        // Change the button's background to the selected color as feedback
        QString style = QString("background-color: %1").arg(color.name());
        colorPickerButton->setStyleSheet(style);

        MRIWidget->colorButton_clicked(color);
    }
}