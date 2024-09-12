#ifndef TRIGGERINGWINDOW_H
#define TRIGGERINGWINDOW_H

#include <QMainWindow>
#include <QMessageBox>
#include <QIntValidator>

#include "../dataHandler/dataHandler.h"

namespace Ui {
class TMSwindow;
}

class TMSwindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit TMSwindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent = nullptr);
    ~TMSwindow();

private slots:
    void on_connectTrigger_clicked();

    void set_enable_UI(bool enable);

    void on_testTrigger_clicked();

    void on_enable_clicked();

    void on_disable_clicked();

    void on_TimeLimitLineEdit_editingFinished();

    void on_setAmplitude_clicked();

    void on_SetMode_clicked();

    void on_RequestInfo_clicked();

private:
    Ui::TMSwindow *ui;

    dataHandler &handler;

    volatile std::sig_atomic_t &signal_received;

    QStringList mode_names ={"Standard", "Power", "Twin", "Dual"};
    QStringList waveform_names ={"StanMonophasicdard", "Biphasic", "Half sine", "Biphasic burst"};
};

#endif // TRIGGERINGWINDOW_H
