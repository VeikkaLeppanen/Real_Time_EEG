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

private slots:
    void updateGraph();

private:

};

#endif // MAINGLWIDGET_H
