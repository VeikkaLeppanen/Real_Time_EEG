#ifndef GLWIDGET_H
#define GLWIDGET_H

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

class Glwidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit Glwidget(QWidget *parent = nullptr);

signals:
    void fetchData();

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

public slots:
    void updateMatrix(const Eigen::MatrixXd &newMatrix, const Eigen::VectorXi &triggers_A, const Eigen::VectorXi &triggers_B);
    void updateChannelDisplayState(std::vector<bool> channelCheckStates);
    void updateGraph();
    void updateChannelNamesQt(QStringList channelNames);
    void updateChannelNamesSTD(std::vector<std::string> channelNames);
    void scaleDrawStateChanged(bool isChecked) { draw_channel_scales = isChecked; }
    void setShowTriggers_A(bool isChecked) { show_triggers_A = isChecked; }
    bool getShowTriggers_A() { return show_triggers_A; }
    void setShowTriggers_B(bool isChecked) { show_triggers_B = isChecked; }
    bool getShowTriggers_B() { return show_triggers_B; }
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

private:
    Eigen::MatrixXd dataMatrix_;
    Eigen::VectorXi triggers_A_;
    Eigen::VectorXi triggers_B_;
    std::vector<bool> channelCheckStates_;
    bool draw_channel_scales = false;
    bool show_triggers_A = true;
    bool show_triggers_B = true;
    bool pause_view = false;
    QStringList channelNames_;

    int matrixCapasity;
    int n_channels;

    double windowLength_seconds;
    int time_line_spacing;
    int totalTimeLines;
    bool drawXaxis = true;
};

#endif // GLWIDGET_H
