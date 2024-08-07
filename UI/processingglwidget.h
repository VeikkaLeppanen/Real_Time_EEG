#ifndef PROCESSINGGLWIDGET_H
#define PROCESSINGGLWIDGET_H

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

class ProcessingGlWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit ProcessingGlWidget(QWidget *parent = nullptr);

signals:
    void fetchData();

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

public slots:
    void updateMatrix(const Eigen::MatrixXd &newMatrix);

    void updateChannelDisplayState(std::vector<bool> channelCheckStates) { channelCheckStates_ = channelCheckStates; }
    void updateGraph();
    void updateChannelNamesQt(QStringList channelNames) { channelNames_ = channelNames; }
    void updateChannelNamesSTD(std::vector<std::string> channelNames);
    void scaleDrawStateChanged(bool isChecked) { draw_channel_scales = isChecked; }

    void setCustomScaleStatus(bool isChecked) { useCustomScale_ = isChecked; }
    void setCustomScaleMin(double min) { min_scale_ = min; }
    void setCustomScaleMax(double max) { max_scale_ = max; }

    bool getCustomScaleStatus() { return useCustomScale_; }
    double getCustomScaleMin() { return min_scale_; }
    double getCustomScaleMax() { return max_scale_; }

private:
    Eigen::MatrixXd dataMatrix_;
    std::vector<bool> channelCheckStates_;
    bool draw_channel_scales = true;
    QStringList channelNames_;

    // Custom scale
    bool useCustomScale_ = false;
    double min_scale_ = -10.0;
    double max_scale_ = 10.0;

    int matrixCapasity_;
    int n_channels_;
};

#endif // PROCESSINGGLWIDGET_H
