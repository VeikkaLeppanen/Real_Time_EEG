#include "triggeringwindow.h"
#include "ui_triggeringwindow.h"

TriggeringWindow::TriggeringWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::TriggeringWindow),
      handler(handler),
      signal_received(signal_received)
{
    ui->setupUi(this);
}

TriggeringWindow::~TriggeringWindow()
{
    delete ui;
}

void TriggeringWindow::on_connectTrigger_clicked()
{
    if(handler.connectTriggerPort()) {
        handler.setTriggerPortStatus(false);
        std::cerr << "Trigger port connection failed" << '\n';
    } else {
        handler.setTriggerPortStatus(true);
        std::cout << "Trigger port connected" << '\n';
    }
}

void TriggeringWindow::on_testTrigger_clicked()
{
    if(handler.getTriggerPortStatus()) {
        handler.trig();
        std::cout << "Trigger sent" << '\n';
    } else {
        std::cerr << "Trigger port not connected" << '\n';
    }
}