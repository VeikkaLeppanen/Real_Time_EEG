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

    ui->lineEditIPI->setValidator(new QIntValidator(1, 3000, this));
    ui->lineEditBAratio->setValidator(new QDoubleValidator(0.2, 5.0, 1, this));

    ui->TimeLimitLineEdit->setText(QString::number(handler.getTriggerTimeLimit()));
    ui->lineEditIPI->setText(QString::number(1));
    ui->lineEditBAratio->setText(QString::number(1));
    ui->checkBoxDelay->setChecked(true);

    ui->ModeLabel->setText("");
    ui->DirectionLabel->setText("");
    ui->WaveformLabel->setText("");
    ui->BurstPulsesLabel->setText("");
    ui->IPILabel->setText("");
    ui->BAratioLabel->setText("");
    ui->EnabledLabel->setText("");
    
    ui->SetModeError->setStyleSheet("QLabel { color : gray; }");
    ui->RequestInfoError->setStyleSheet("QLabel { color : gray; }");

    setWindowTitle("TMS Window");
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






void TriggeringWindow::on_SetMode_clicked()
{
    if(handler.getTriggerConnectStatus()) {
        ui->SetModeError->setText("Set mode requested");

        // Get values from comboboxes
        int mode = ui->comboBoxMode->currentIndex();
        int direction = ui->comboBoxDirection->currentIndex();
        int waveform = ui->comboBoxWaveform->currentIndex();
        int burst_pulses = ui->comboBoxBurstPulses->currentText().toInt();

        // Get values from line edits
        int ipi = ui->lineEditIPI->text().toInt();
        double ba_ratio = ui->lineEditBAratio->text().toDouble();

        // Get the state of the checkbox
        bool delay = ui->checkBoxDelay->isChecked();

        ui->SetModeError->setText("Sending set mode request");
        handler.magPro_set_mode(mode, direction, waveform, burst_pulses, ipi, ba_ratio, delay);

        ui->SetModeError->setText("magPro mode set");
    } else {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected.");
    }
}


void TriggeringWindow::on_RequestInfo_clicked()
{
    if(handler.getTriggerConnectStatus()) {
        ui->RequestInfoError->setText("Mode info requested");

        handler.magPro_request_mode_info();
        ui->RequestInfoError->setText("Mode info request sent");

        int mode;
        int direction;
        int waveform;
        int burst_pulses;
        float ipi;
        float ba_ratio;
        bool enabled;
        handler.get_mode_info(mode, direction, waveform, burst_pulses, ipi, ba_ratio, enabled);
        ui->RequestInfoError->setText("Mode info retrieved");

        ui->ModeLabel->setText(mode_names.at(mode));
        ui->DirectionLabel->setText(direction ? "Reverse" : "Normal");
        ui->WaveformLabel->setText(waveform_names.at(waveform));
        ui->BurstPulsesLabel->setText(QString::number(burst_pulses));
        ui->IPILabel->setText(QString::number(ipi));
        ui->BAratioLabel->setText(QString::number(ba_ratio));
        ui->EnabledLabel->setText(enabled ? "Enabled" : "Disabled");

    } else {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected.");
    }
}

