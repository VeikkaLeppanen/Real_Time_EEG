#ifndef MAINGLWIDGET_H
#define MAINGLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QTimer>
#include <QLabel>
#include <algorithm>
#include "../dataHandler/dataHandler.h"
#include <QPainter>
#include <QFont>
#include <QString>

#include <algorithm>

class MainGlWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit MainGlWidget(QWidget *parent = nullptr);

signals:
    void fetchData();

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

public slots:
    void updateMatrix(Eigen::MatrixXd &newMatrix) {
    // Check if dimensions differ
    if (dataMatrix_.rows() != newMatrix.rows() || dataMatrix_.cols() != newMatrix.cols()) {
        // Resize dataMatrix_ to match the dimensions of newMatrix
        dataMatrix_.resize(newMatrix.rows(), newMatrix.cols());
    }

    // Assign the new matrix
    dataMatrix_ = newMatrix;

    // Update other attributes based on the new matrix
    matrixCapasity_ = newMatrix.cols(); 
    n_channels_ = newMatrix.rows();
}

    void updateChannelDisplayState(std::vector<bool> channelCheckStates) { channelCheckStates_ = channelCheckStates; }
    void updateGraph();
    void updateChannelNamesQt(QStringList channelNames) { channelNames_ = channelNames; }
    void updateChannelNamesSTD(std::vector<std::string> channelNames) {
        QStringList newNames;
        for(size_t i = 0; i < channelNames.size(); i++) {
            newNames.append(QString::fromStdString(channelNames[i])); 
        }
        channelNames_ = newNames;
    }
    void scaleDrawStateChanged(bool isChecked) { draw_channel_scales = isChecked; }

private:
    Eigen::MatrixXd dataMatrix_;
    std::vector<bool> channelCheckStates_;
    bool draw_channel_scales = false;
    QStringList channelNames_;

    int matrixCapasity_;
    int n_channels_;
};

#endif // MAINGLWIDGET_H
