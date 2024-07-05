#include "triggeringwindow.h"
#include "ui_triggeringwindow.h"

TriggeringWindow::TriggeringWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::TriggeringWindow),
      handler(handler),
      signal_received(signal_received)
{
    ui->setupUi(this);

    ui->amplitudeLineEdit->setValidator(new QIntValidator(0, 100, this));
    ui->TimeLimitLineEdit->setValidator(new QIntValidator(100, 100000, this));

    ui->TimeLimitLineEdit->setText(QString::number(handler.getTriggerTimeLimit()));
}

TriggeringWindow::~TriggeringWindow()
{
    delete ui;
}

void TriggeringWindow::on_connectTrigger_clicked()
{
    if(handler.connectTriggerPort()) {
        handler.setTriggerConnectStatus(false);
        std::cerr << "Trigger port connection failed" << '\n';
    } else {
        handler.setTriggerConnectStatus(true);
        std::cout << "Trigger port connected" << '\n';
    }
}

void TriggeringWindow::on_testTrigger_clicked()
{
    if(handler.getTriggerEnableStatus()) {
        handler.send_trigger();
        std::cout << "Trigger sent" << '\n';
    } else {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not enabled.");
    }
}
void TriggeringWindow::on_enable_clicked()
{
    if(handler.getTriggerConnectStatus()) {
        handler.set_enable(true);
        std::cout << "Trigger port enabled" << '\n';
    } else {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected.");
    }
}


void TriggeringWindow::on_disable_clicked()
{
    if(handler.getTriggerConnectStatus()) {
        handler.set_enable(false);
        std::cout << "Trigger port disabled" << '\n';
    } else {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected.");
    }
}


void TriggeringWindow::on_TimeLimitLineEdit_editingFinished()
{
    bool ok;
    int value = ui->TimeLimitLineEdit->text().toInt(&ok);
    if (ok) {
        handler.setTriggerTimeLimit(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 100 and 100000.");
    }
}


void TriggeringWindow::on_setAmplitude_clicked()
{
    bool ok;
    int value = ui->amplitudeLineEdit->text().toInt(&ok);
    if (!ok) {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 100.");
    } else if (!handler.getTriggerEnableStatus()) {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not enabled.");
    } else {
        handler.set_amplitude(value);
    }
}



