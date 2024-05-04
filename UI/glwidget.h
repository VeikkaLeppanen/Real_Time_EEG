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
    void updateMatrix(Eigen::MatrixXd newMatrix) { newMatrix_ = newMatrix; }
    void updateGraph();
    void updateChannelNames(std::vector<std::string> channelNames) {
        QStringList newNames;
        for(size_t i = 0; i < channelNames.size(); i++) {
            newNames.append(QString::fromStdString(channelNames[i])); 
        }
        channelNames_ = newNames;
    }

private:
    Eigen::MatrixXd newMatrix_;  // Use Eigen::VectorXf if you need a vector of floats
    QStringList channelNames_;

};

#endif // GLWIDGET_H
