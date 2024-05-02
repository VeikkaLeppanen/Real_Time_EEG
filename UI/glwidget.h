#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QTimer>
#include <algorithm>
#include "../dataHandler/dataHandler.h"

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
    void updateMatrix(Eigen::MatrixXd newMatrix) { dataMatrix = newMatrix; }
    void updateGraph();

private:
    Eigen::MatrixXd dataMatrix;  // Use Eigen::VectorXf if you need a vector of floats

};

#endif // GLWIDGET_H
