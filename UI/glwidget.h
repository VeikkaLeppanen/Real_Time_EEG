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

private:
    Eigen::MatrixXd dataMatrix_;
    Eigen::VectorXi triggers_A_;
    Eigen::VectorXi triggers_B_;
    std::vector<bool> channelCheckStates_;
    bool draw_channel_scales = false;
    QStringList channelNames_;

    int matrixCapasity;
    int n_channels;
};

#endif // GLWIDGET_H
