#include "TMSwindow.h"
#include "../ui_TMSwindow.h"

TMSwindow::TMSwindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent)
    : QMainWindow(parent),
      ui(new Ui::TMSwindow),
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

    ui->ModeLabel->setText("");
    ui->DirectionLabel->setText("");
    ui->WaveformLabel->setText("");
    ui->BurstPulsesLabel->setText("");
    ui->IPILabel->setText("");
    ui->BAratioLabel->setText("");
    ui->EnabledLabel->setText("");
    
    ui->SetModeError->setStyleSheet("QLabel { color : gray; }");
    ui->RequestInfoError->setStyleSheet("QLabel { color : gray; }");

    if (handler.getTriggerConnectStatus()) set_enable_UI(true);
    else set_enable_UI(false);

    setWindowTitle("TMS Window");
}

TMSwindow::~TMSwindow()
{
    delete ui;
}

void TMSwindow::on_connectTrigger_clicked()
{
    switch (connectionType) {
        case COM:
            if(handler.connectTriggerPort()) {
                handler.setTriggerConnectStatus(false);
                set_enable_UI(false);
                std::cerr << "Trigger port connection failed" << '\n';
            } else {
                handler.setTriggerConnectStatus(true);
                set_enable_UI(true);
                on_RequestInfo_clicked();
                handler.set_amplitude(0);
                std::cout << "Trigger port connected" << '\n';
            }

            break;

        case TTL:
            if(handler.connectTriggerPort_TTL()) {
                handler.setTriggerConnectStatus(false);
                set_enable_UI(false);
                std::cerr << "Trigger port connection failed" << '\n';
            } else {
                handler.setTriggerConnectStatus(true);
                set_enable_UI(true);
                std::cout << "Trigger port connected" << '\n';
            }

            break;
    }
}

void TMSwindow::set_enable_UI(bool enable) 
{
    switch (connectionType) {
        case COM:
            ui->comboBoxBurstPulses->setEnabled(enable);
            ui->comboBoxDirection->setEnabled(enable);
            ui->comboBoxMode->setEnabled(enable);
            ui->comboBoxWaveform->setEnabled(enable);

            ui->enable->setEnabled(enable);
            ui->disable->setEnabled(enable);
            ui->setAmplitude->setEnabled(enable);
            ui->testTrigger->setEnabled(enable);
            ui->SetMode->setEnabled(enable);
            ui->RequestInfo->setEnabled(enable);

            ui->amplitudeLineEdit->setEnabled(enable);
            ui->lineEditIPI->setEnabled(enable);
            ui->lineEditBAratio->setEnabled(enable);
            ui->TimeLimitLineEdit->setEnabled(enable);

            ui->ModeLabel->setEnabled(enable);
            ui->DirectionLabel->setEnabled(enable);
            ui->WaveformLabel->setEnabled(enable);
            ui->BurstPulsesLabel->setEnabled(enable);
            ui->IPILabel->setEnabled(enable);
            ui->BAratioLabel->setEnabled(enable);
            ui->EnabledLabel->setEnabled(enable);

            ui->SetModeError->setEnabled(enable);
            ui->RequestInfoError->setEnabled(enable);

            break;

        case TTL:
            ui->comboBoxBurstPulses->setEnabled(false);
            ui->comboBoxDirection->setEnabled(false);
            ui->comboBoxMode->setEnabled(false);
            ui->comboBoxWaveform->setEnabled(false);

            ui->enable->setEnabled(enable);
            ui->disable->setEnabled(enable);
            ui->setAmplitude->setEnabled(false);
            ui->testTrigger->setEnabled(enable);
            ui->SetMode->setEnabled(false);
            ui->RequestInfo->setEnabled(false);

            ui->amplitudeLineEdit->setEnabled(false);
            ui->lineEditIPI->setEnabled(false);
            ui->lineEditBAratio->setEnabled(false);
            ui->TimeLimitLineEdit->setEnabled(enable);

            ui->ModeLabel->setEnabled(false);
            ui->DirectionLabel->setEnabled(false);
            ui->WaveformLabel->setEnabled(false);
            ui->BurstPulsesLabel->setEnabled(false);
            ui->IPILabel->setEnabled(false);
            ui->BAratioLabel->setEnabled(false);
            ui->EnabledLabel->setEnabled(false);

            ui->SetModeError->setEnabled(false);
            ui->RequestInfoError->setEnabled(false);

            break;
    }
}

void TMSwindow::on_testTrigger_clicked()
{
    switch (connectionType) {
        case COM:
            if(handler.getTriggerEnableStatus()) {
                handler.send_trigger();
                std::cout << "Trigger sent" << '\n';
            } else { QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected."); }

            break;

        case TTL:
            if(handler.getTriggerEnableStatus()) {
                handler.send_trigger_TTL();
                std::cout << "Trigger sent" << '\n';
            } else { QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected."); }

            break;
    }
}
void TMSwindow::on_enable_clicked()
{
    switch (connectionType) {
        case COM:
            if(handler.getTriggerConnectStatus()) {
                handler.set_enable(true);
                std::cout << "Trigger port enabled" << '\n';
            } else { QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected."); }

            break;

        case TTL:
            if(handler.getTriggerConnectStatus()) {
                handler.set_enable_TTL(true);
                std::cout << "Trigger port enabled" << '\n';
            } else { QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected."); }

            break;
    }
}


void TMSwindow::on_disable_clicked()
{
    switch (connectionType) {
        case COM:
            if(handler.getTriggerConnectStatus()) {
                handler.set_enable(false);
                std::cout << "Trigger port disabled" << '\n';
            } else { QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected."); }

            break;

        case TTL:
            if(handler.getTriggerConnectStatus()) {
                handler.set_enable_TTL(false);
                std::cout << "Trigger port disabled" << '\n';
            } else { QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected."); }

            break;
    }
}


void TMSwindow::on_TimeLimitLineEdit_editingFinished()
{
    bool ok;
    int value = ui->TimeLimitLineEdit->text().toInt(&ok);
    if (ok && value >= 100 && value <= 100000) {
        handler.setTriggerTimeLimit(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 100 and 100000.");
    }
}


void TMSwindow::on_setAmplitude_clicked()
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






void TMSwindow::on_SetMode_clicked()
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
        bool delay = true;

        ui->SetModeError->setText("Sending set mode request");
        handler.magPro_set_mode(mode, direction, waveform, burst_pulses, ipi, ba_ratio, delay);

        ui->SetModeError->setText("magPro mode set");
    } else {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected.");
    }
}


void TMSwindow::on_RequestInfo_clicked()
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

        // Set the values in the labels
        ui->ModeLabel->setText(mode_names.at(mode));
        ui->DirectionLabel->setText(direction ? "Reverse" : "Normal");
        ui->WaveformLabel->setText(waveform_names.at(waveform));
        ui->BurstPulsesLabel->setText(QString::number(burst_pulses));
        ui->IPILabel->setText(QString::number(ipi));
        ui->BAratioLabel->setText(QString::number(ba_ratio));
        ui->EnabledLabel->setText(enabled ? "Enabled" : "Disabled");

        // Set the values in the comboboxes
        ui->comboBoxMode->setCurrentIndex(mode);
        ui->comboBoxDirection->setCurrentIndex(direction);
        ui->comboBoxWaveform->setCurrentIndex(waveform);
        ui->comboBoxBurstPulses->setCurrentIndex(burst_pulses);

        // Set the values in the line edits
        ui->lineEditIPI->setText(QString::number(ipi, 'f', 2));
        ui->lineEditBAratio->setText(QString::number(ba_ratio, 'f', 2));

    } else {
        QMessageBox::warning(this, "Trigger Port Error", "Trigger port not connected.");
    }
}


void TMSwindow::on_comboBox_connectionType_currentIndexChanged(int index)
{
    if (index == 0) { connectionType = COM; }
    else if (index == 1) { connectionType = TTL; }
    handler.setTMSConnectionType(connectionType);
}
