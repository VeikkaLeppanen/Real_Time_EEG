#ifndef PROCESSINGGLWIDGET_H
#define PROCESSINGGLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QTimer>
#include <QLabel>
#include <algorithm>
#include <QPainter>
#include <QFont>
#include <QString>
#include <algorithm>
#include <QMouseEvent>  // Add this line at the top of your processingglwidget.cpp

#include "../dataHandler/dataHandler.h"
#include "customTooltip.h"

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
                      const Eigen::VectorXi &triggers_out, 
                      const Eigen::VectorXd &time_stamps,
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
    void setShowTriggers_out(bool isChecked) { show_triggers_out = isChecked; }
    bool getShowTriggers_out() { return show_triggers_out; }
    void switchPause() { pause_view = !pause_view; }

    void updateWindowLength_seconds(double windowLength) {
        windowLength_seconds = windowLength;
        totalTimeLines = (windowLength_seconds * 1000) / time_line_spacing;
    }

    void updateTLineSpacing(int LineSpacing) {
        time_line_spacing = LineSpacing;
        totalTimeLines = (windowLength_seconds * 1000) / time_line_spacing;
    }
    void setDrawXaxis(bool isChecked) { drawXaxis = isChecked; }
    bool getDrawXaxis() { return drawXaxis; }
    void mouseMoveEvent(QMouseEvent *event);
    QString getDataValue(const QPoint &pos);
    void leaveEvent(QEvent *event);
    
private:
    std::mutex dataMutex;
    
    Eigen::MatrixXd dataMatrix_;
    Eigen::VectorXi triggers_A_;
    Eigen::VectorXi triggers_B_;
    Eigen::VectorXi triggers_out_;
    Eigen::VectorXd time_stamps_;
    std::vector<bool> channelCheckStates_;
    bool draw_channel_scales = true;
    bool show_triggers_A = true;
    bool show_triggers_B = true;
    bool show_triggers_out = true;
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

    double windowLength_seconds;
    int time_line_spacing;
    int totalTimeLines;
    bool drawXaxis = true;

    CustomTooltip *tooltipWidget = nullptr;
    QPoint currentMousePosition;
    int visibleSampleCount;
    float currentTimePosition;
};

#endif // PROCESSINGGLWIDGET_H
