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
    void updateMatrix(Eigen::MatrixXd newMatrix) { 
        newMatrix_ = newMatrix;
        matrixCapasity = newMatrix.cols(); 
        n_channels = newMatrix.rows();
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

private:
    Eigen::MatrixXd newMatrix_;
    std::vector<bool> channelCheckStates_;
    QStringList channelNames_;

    int matrixCapasity;
    int n_channels;
};

#endif // GLWIDGET_H
