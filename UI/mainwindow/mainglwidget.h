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

#include "../dataHandler/dataHandler.h"
#include <image/image.h>
#include <image/image_operators.h>
#include <image/orientation.h>
#include <Eigen/Dense>

class MainGlWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    enum class SignalType {
        EEG,
        MRI
    };

    explicit MainGlWidget(QWidget *parent = nullptr);
    void setDataHandler(dataHandler* handler) { dataHandler_ = handler; }

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

signals:

public slots:
    void setSignalSource(int sourceIndex);
    void updateData(const Eigen::VectorXd &newData, 
                    const Eigen::VectorXi &triggers_A, 
                    const Eigen::VectorXi &triggers_B, 
                    const Eigen::VectorXi &triggers_out, 
                    const Eigen::VectorXd &time_stamps,
                             const size_t &data_index,
                              std::string source_name);
    void setTriggerAVisible(bool visible) { show_triggers_A = visible; }
    void setTriggerBVisible(bool visible) { show_triggers_B = visible; }
    void setTriggerOutVisible(bool visible) { show_trigger_out = visible; }
    void setViewerType(int type) { viewerType = type; update(); }

private slots:
    void updateGraph();

private:
    std::mutex dataMutex;
    
    Eigen::VectorXd dataVector_;
    Eigen::VectorXi triggers_A_;
    Eigen::VectorXi triggers_B_;
    Eigen::VectorXi triggers_out_;
    Eigen::VectorXd time_stamps_;
    int data_index_;
    bool draw_channel_scales = false;
    bool show_triggers_A = false;
    bool show_triggers_B = false;
    bool show_trigger_out = false;
    bool pause_view = false;
    QString source_name_;

    int vectorCapasity;

    double windowLength_seconds;
    int time_line_spacing;
    int totalTimeLines;
    bool drawXaxis = true;

    enum class ScalingMethod {
        MINMAX_COEFF,
        MEAN_STD
    };
    
    std::pair<double, double> calculateDataRange(const Eigen::VectorXd& data, ScalingMethod method = ScalingMethod::MINMAX_COEFF);

    void disconnectCurrentSource();
    void connectToSource(const QString &source);
    
    dataHandler* dataHandler_ = nullptr;
    QMetaObject::Connection currentConnection;
    QString currentSource;

    int viewerType = 0;  // 0 for Signal viewer, 1 for Response monitor

    void paintSignalViewer();  // New function for signal viewer mode
    void paintResponseMonitor();  // New function for response monitor mode

    int trigger_position_tolerance = 50;
    int response_window_size = 10000;  // Fixed window size for response monitor
    float samplingRate = 5000.0f;  // Hz (adjust as needed)

    SignalType currentSignalType = SignalType::EEG;  // Default to EEG
};

#endif // MAINGLWIDGET_H
