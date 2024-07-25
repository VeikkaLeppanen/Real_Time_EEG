#ifndef TRIGGERINGWINDOW_H
#define TRIGGERINGWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <QIntValidator>

#include "../dataHandler/dataHandler.h"

namespace Ui {
class TriggeringWindow;
}

class TriggeringWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit TriggeringWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent = nullptr);
    ~TriggeringWindow();

private slots:
    void on_connectTrigger_clicked();

    void on_testTrigger_clicked();

    void on_enable_clicked();

    void on_disable_clicked();

    void on_TimeLimitLineEdit_editingFinished();

    void on_setAmplitude_clicked();

    void on_SetMode_clicked();

    void on_RequestInfo_clicked();

private:
    Ui::TriggeringWindow *ui;

    dataHandler &handler;

    volatile std::sig_atomic_t &signal_received;

    QStringList mode_names ={"Standard", "Power", "Twin", "Dual"};
    QStringList waveform_names ={"StanMonophasicdard", "Biphasic", "Half sine", "Biphasic burst"};
};

#endif // TRIGGERINGWINDOW_H
