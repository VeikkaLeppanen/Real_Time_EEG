#ifndef MRIWINDOW_H
#define MRIWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include <QToolBar>
#include <QStatusBar>
#include <QMenuBar>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QListWidget>
#include <QListWidgetItem>
#include <QSlider>
#include <QColorDialog>
#include <QComboBox>
#include <QButtonGroup>
#include <QPushButton>
#include <QSpinBox>
#include <QLabel>
#include <QMessageBox>
#include <QFileDialog>
#include <QLineEdit>
#include "openglMRIwidget.h"  // Include your custom widget header
#include "../dataHandler/dataHandler.h"

class mriWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit mriWindow(dataHandler &handler, QWidget *parent = nullptr);
    ~mriWindow();

    OpenGLMRIWidget* getGlWidget() { return MRIWidget; }

public slots:

    void ROI_update(const std::vector<std::string> &names, std::vector<bool> ROI_toggleStatus) {
        // Clear the existing items in the list
        toggleList->clear();
        
        // Add new items from the names vector
        for (int i = 0; i < names.size(); ++i) {
            QListWidgetItem *item = new QListWidgetItem(QString::fromStdString(names[i]), toggleList);
            item->setFlags(item->flags() | Qt::ItemIsUserCheckable); // Make item checkable
            item->setCheckState(ROI_toggleStatus[i] ? Qt::Checked : Qt::Unchecked); // Set the initial state to unchecked
            toggleList->addItem(item);
        }
        if (toggleList->count() <= 0) {
            editButton->setStyleSheet("");
            editButton->setChecked(false);
        }
    }

private slots:
    
    void MRI_T1_load_clicked();
    void fMRI_load_clicked();
    void setWatchFolder();
    
    void toggleToolbar() {
        if (mainToolbar)
            mainToolbar->setVisible(!mainToolbar->isVisible());
    }
    void updateTargetROI(QListWidgetItem *item) { MRIWidget->updateTargetROI(toggleList->row(item)); };

    void editButton_toggled(bool checked) {
        if (toggleList->count() > 0) {
            if (checked) editButton->setStyleSheet("background-color: darkgreen;");
            else editButton->setStyleSheet("");

            MRIWidget->editButton_toggled(checked);
        } else {
            editButton->setStyleSheet("");
            editButton->setChecked(false);
            MRIWidget->editButton_toggled(false);
            QMessageBox::warning(this, "Warning", "No ROI loaded.");
        }
    }

    void updateToggleStatus();
    void loadButton_clicked();
    void colorButton_clicked();

private:
    dataHandler &handler;  // Add handler reference
    OpenGLMRIWidget* MRIWidget;  // Assuming you have a Glwidget class
    QWidget* centralWidget;
    QVBoxLayout* mainLayout;
    QToolBar *mainToolbar;

    QPushButton *editButton;
    QListWidget *toggleList;
    QPushButton *colorPickerButton;
};

#endif // MRIWINDOW_H
