#ifndef MAINGLWIDGET_H
#define MAINGLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QTimer>
#include <iostream>
#include <float.h>
#include <unordered_set>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QPainter>
#include <QFont>
#include <image/image.h>
#include <image/image_operators.h>
#include <image/orientation.h>
#include <Eigen/Dense>

class MainGlWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit MainGlWidget(QWidget *parent = nullptr);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

signals:

public slots:
    void setSignalSource(const QString &source);
    void updateData(const Eigen::VectorXi &newData, 
                    const Eigen::VectorXi &triggers_A, 
                    const Eigen::VectorXi &triggers_B, 
                    const Eigen::VectorXd &time_stamps,
                              std::string source_name);

private slots:
    void updateGraph();

private:
    std::mutex dataMutex;
    
    Eigen::VectorXi dataVector_;
    Eigen::VectorXi triggers_A_;
    Eigen::VectorXi triggers_B_;
    Eigen::VectorXd time_stamps_;
    bool draw_channel_scales = false;
    bool show_triggers_A = true;
    bool show_triggers_B = true;
    bool pause_view = false;
    QString source_name_;

    int vectorCapasity;

    double windowLength_seconds;
    int time_line_spacing;
    int totalTimeLines;
    bool drawXaxis = true;

};

#endif // MAINGLWIDGET_H
