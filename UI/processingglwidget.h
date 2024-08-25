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
    void updateMatrix(const Eigen::MatrixXd &newMatrix, 
                      const Eigen::VectorXi &triggers_A, 
                      const Eigen::VectorXi &triggers_B, 
                                        int numPastElements, 
                                        int numFutureElements);

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
    void setShowTriggers_A(bool isChecked) { show_triggers_A = isChecked; }
    bool getShowTriggers_A() { return show_triggers_A; }
    void setShowTriggers_B(bool isChecked) { show_triggers_B = isChecked; }
    bool getShowTriggers_B() { return show_triggers_B; }
    void switchPause() { pause_view = !pause_view; }

private:
    Eigen::MatrixXd dataMatrix_;
    Eigen::VectorXi triggers_A_;
    Eigen::VectorXi triggers_B_;
    std::vector<bool> channelCheckStates_;
    bool draw_channel_scales = true;
    bool show_triggers_A = true;
    bool show_triggers_B = true;
    bool pause_view = false;
    QStringList channelNames_;

    int numPastElements_;
    int numFutureElements_;

    // Custom scale
    bool useCustomScale_ = false;
    double min_scale_ = -10.0;
    double max_scale_ = 10.0;

    int matrixCapasity_;
    int n_channels_;
};

#endif // PROCESSINGGLWIDGET_H
