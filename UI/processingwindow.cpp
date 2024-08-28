#include "processingwindow.h"
#include "ui_processingwindow.h"
#include <QThread>

ProcessingWindow::ProcessingWindow(dataHandler &handler, 
                    volatile std::sig_atomic_t &processingWorkerRunning, 
                               Eigen::MatrixXd &processed_data, 
                                       QWidget *parent)
    : QMainWindow(parent), 
      ui(new Ui::ProcessingWindow),
      handler(handler),
      processingWorkerRunning(processingWorkerRunning),
      processed_data(processed_data)
{
    ui->setupUi(this);

    // Initialize values for lineEdits
    ui->edge->setText(QString::number(phaseEstParams.edge));
    ui->modelOrder->setText(QString::number(phaseEstParams.modelOrder));
    ui->hilbertLength->setText(QString::number(phaseEstParams.hilbertWinLength));
    ui->stimulationTarget->setText(QString::number(phaseEstParams.stimulation_target));
    ui->phaseShift->setText(QString::number(phaseEstParams.phase_shift));
    ui->lineEdit_XaxisSpacing->setText(QString::number(500));

    ui->checkBox_triggers_A->setStyleSheet("QCheckBox { color : blue; }");
    ui->checkBox_triggers_B->setStyleSheet("QCheckBox { color : green; }");
    ui->checkBox_triggers_out->setStyleSheet("QCheckBox { color : cyan; }");
    setWindowTitle("Phase Estimation Window");
    resize(1280, 720);

    processingglWidget = ui->processingGlWidget;
    if (processingglWidget) {
        connect(this, &ProcessingWindow::setCustomScaleStatus, processingglWidget, &ProcessingGlWidget::setCustomScaleStatus);
        connect(this, &ProcessingWindow::setCustomScaleMin, processingglWidget, &ProcessingGlWidget::setCustomScaleMin);
        connect(this, &ProcessingWindow::setCustomScaleMax, processingglWidget, &ProcessingGlWidget::setCustomScaleMax);

        connect(this, &ProcessingWindow::getCustomScaleStatus, processingglWidget, &ProcessingGlWidget::getCustomScaleStatus);
        connect(this, &ProcessingWindow::getCustomScaleMin, processingglWidget, &ProcessingGlWidget::getCustomScaleMin);
        connect(this, &ProcessingWindow::getCustomScaleMax, processingglWidget, &ProcessingGlWidget::getCustomScaleMax);
        connect(this, &ProcessingWindow::setShowTriggers_A, processingglWidget, &ProcessingGlWidget::setShowTriggers_A);
        connect(this, &ProcessingWindow::setShowTriggers_B, processingglWidget, &ProcessingGlWidget::setShowTriggers_B);
        connect(this, &ProcessingWindow::setShowTriggers_out, processingglWidget, &ProcessingGlWidget::setShowTriggers_out);
        connect(this, &ProcessingWindow::switchPause, processingglWidget, &ProcessingGlWidget::switchPause);

        connect(this, &ProcessingWindow::updateWidgetChannelNames, processingglWidget, &ProcessingGlWidget::updateChannelNamesSTD);
        connect(this, &ProcessingWindow::setDrawXaxis, processingglWidget, &ProcessingGlWidget::setDrawXaxis);
        connect(this, &ProcessingWindow::updateTLineSpacing, processingglWidget, &ProcessingGlWidget::updateTLineSpacing);

        ui->checkBox->setChecked(emit getCustomScaleStatus());
        ui->scaleMin->setText(QString::number(emit getCustomScaleMin()));
        ui->scaleMax->setText(QString::number(emit getCustomScaleMax()));
        ui->checkBox_triggers_A->setChecked(processingglWidget->getShowTriggers_A());
        ui->checkBox_triggers_B->setChecked(processingglWidget->getShowTriggers_B());
        ui->checkBox_triggers_out->setChecked(processingglWidget->getShowTriggers_out());
        ui->checkBox_Xaxis->setChecked(processingglWidget->getDrawXaxis());
    } else {
        // Error handling if glWidget is not found
        qWarning("Glwidget not found in UI!");
    }
}

void ProcessingWindow::updateSpatialChannelNames(std::vector<std::string> names) 
{
    spatial_channel_names = names;
}

ProcessingWindow::~ProcessingWindow()
{
    delete ui;
}

void ProcessingWindow::on_edge_editingFinished()
{
    bool ok;
    int value = ui->edge->text().toInt(&ok);
    if (ok) {
        phaseEstParams.edge = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_modelOrder_editingFinished()
{
    bool ok;
    int value = ui->modelOrder->text().toInt(&ok);
    if (ok) {
        phaseEstParams.modelOrder = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_hilbertLength_editingFinished()
{
    bool ok;
    int value = ui->hilbertLength->text().toInt(&ok);
    if (ok) {
        phaseEstParams.hilbertWinLength = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_stimulationTarget_editingFinished()
{
    bool ok;
    double value = ui->stimulationTarget->text().toDouble(&ok);
    if (ok) {
        phaseEstParams.stimulation_target = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_phaseShift_editingFinished()
{
    bool ok;
    int value = ui->phaseShift->text().toInt(&ok);
    if (ok) {
        phaseEstParams.phase_shift = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}




void ProcessingWindow::on_scaleMax_editingFinished()
{
    bool ok;
    double value = ui->scaleMax->text().toDouble(&ok);
    if (ok) {
        emit setCustomScaleMax(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_scaleMin_editingFinished()
{
    bool ok;
    double value = ui->scaleMin->text().toDouble(&ok);
    if (ok) {
        emit setCustomScaleMin(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_checkBox_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setCustomScaleStatus(isChecked);
}


void ProcessingWindow::on_Filter_checkbox_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setFilterState(isChecked);
}


void ProcessingWindow::on_checkBox_Channels_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setEEGViewState(isChecked);
}


void ProcessingWindow::on_checkBox_PhaseTargeting_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setPhaseTargetingState(isChecked);
}

void ProcessingWindow::on_checkBox_Stimulation_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    handler.setTriggerEnableStatus(isChecked);
}

void ProcessingWindow::on_checkBox_phaseEstimate_stateChanged(int arg1)
{
    if (handler.channelNamesSet()) {
        sampling_rate_ = handler.getSamplingRate();
        bool isChecked = (arg1 == Qt::Checked);
        emit setphaseEstimateState(isChecked);
        
        ui->comboBox_spatialTarget->clear(); // Clear the combo box before updating it

        for (const auto &name : spatial_channel_names)
        {
            ui->comboBox_spatialTarget->addItem(QString::fromStdString(name)); // Convert std::string to QString and add to combo box
        }

        setupComboBox_OuterElectrode(0);
    
    } else {
        QMessageBox::warning(this, "Channel Naming Error", "Please load channel map in the EEG window.");
    }
}


void ProcessingWindow::on_setParamsButton_clicked()
{
    emit setPhaseEstParams(phaseEstParams);
}

void ProcessingWindow::on_comboBox_spatialTarget_currentIndexChanged(int index)
{
    emit setSpatilaTargetChannel(index);
    setupComboBox_OuterElectrode(index);
}

void ProcessingWindow::setupComboBox_OuterElectrode(int spatial_target_index) {

    // Setup outer channel names
    outer_channel_names.clear();
    bool skip_first = true;

    for (int i = 0; i < spatial_channel_names.size(); i++)
    {
        if (i != spatial_target_index) {
            outer_channel_names.append(QString::fromStdString(spatial_channel_names[i]));
        }
    }

    // Setup combobox
    QStandardItemModel* model = new QStandardItemModel(numOuterElectrodes, 1, this);
    outerElectrodeCheckStates_.resize(numOuterElectrodes, true);

    for (int i = 0; i < numOuterElectrodes; ++i) {
        QStandardItem* item = new QStandardItem(outer_channel_names.at(i));
        item->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
        item->setData(Qt::Checked, Qt::CheckStateRole);
        model->setItem(i, 0, item);
    }

    ui->comboBox_outElectrodes->setModel(model);
    connect(model, &QStandardItemModel::itemChanged, this, &ProcessingWindow::handleCheckboxChange);
}

void ProcessingWindow::handleCheckboxChange(QStandardItem* item) {
    int row = item->row();
    bool isChecked = item->checkState() == Qt::Checked;
    outerElectrodeCheckStates_[row] = isChecked;
    
    emit outerElectrodesStateChanged(outerElectrodeCheckStates_);
}

void ProcessingWindow::on_pushButton_pause_view_clicked()
{
    emit switchPause();
}

void ProcessingWindow::on_checkBox_triggers_A_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_A(isChecked);
}

void ProcessingWindow::on_checkBox_triggers_B_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_B(isChecked);
}

void ProcessingWindow::on_checkBox_phaseError_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setPhaseError(isChecked);
}

void ProcessingWindow::on_checkBox_Xaxis_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setDrawXaxis(isChecked);
}

void ProcessingWindow::on_lineEdit_XaxisSpacing_editingFinished()
{
    bool ok;
    int value = ui->lineEdit_XaxisSpacing->text().toInt(&ok);
    if (ok) {
        emit updateTLineSpacing(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void ProcessingWindow::on_checkBox_triggers_out_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_out(isChecked);
}


void ProcessingWindow::on_comboBox_phaseError_currentIndexChanged(int index)
{
    emit setPhaseErrorType(index);
}

