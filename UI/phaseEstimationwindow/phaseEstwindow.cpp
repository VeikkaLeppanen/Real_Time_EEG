#include "phaseEstwindow.h"
#include "../ui_phaseEstwindow.h"
#include <QThread>
// #include "phaseEstwindow.moc"

phaseEstwindow::phaseEstwindow(dataHandler &handler, 
                    volatile std::sig_atomic_t &processingWorkerRunning, 
                               Eigen::MatrixXd &processed_data, 
                                       QWidget *parent)
    : QMainWindow(parent), 
      ui(new Ui::phaseEstwindow),
      handler(handler),
      processingWorkerRunning(processingWorkerRunning),
      processed_data(processed_data)
{
    ui->setupUi(this);

    // Initialize values for lineEdits
    ui->SNRmax->setText(QString::number(0.0f));
    ui->SNRmax->setAttribute(Qt::WA_Hover);         // Enable hover events
    ui->SNRmax->setMouseTracking(true);             // Enable mouse tracking
    ui->SNRmax->setAttribute(Qt::WA_AlwaysShowToolTips);  // Always show tooltips, even if the widget is disabled

    ui->SNRtreshold->setText(QString::number(phaseEstParams.SNR_threshold));
    ui->edge->setText(QString::number(phaseEstParams.edge));
    ui->modelOrder->setText(QString::number(phaseEstParams.modelOrder));
    ui->hilbertLength->setText(QString::number(phaseEstParams.hilbertWinLength));
    ui->stimulationTarget->setText(QString::number(phaseEstParams.stimulation_target));
    ui->phaseShift->setText(QString::number(phaseEstParams.phase_shift));
    ui->lineEdit_XaxisSpacing->setText(QString::number(500));
    ui->checkBox_Stimulation->setChecked(handler.getTriggerEnableStatus());

    ui->tabWidget->setTabText(0, "Phase Estimation");
    ui->tabWidget->setTabText(1, "View");
    ui->checkBox_triggers_A->setStyleSheet("QCheckBox { color : blue; }");
    ui->checkBox_triggers_B->setStyleSheet("QCheckBox { color : green; }");
    ui->checkBox_triggers_out->setStyleSheet("QCheckBox { color : cyan; }");
    setWindowTitle("Phase Estimation Window");
    resize(1280, 720);

    processingglWidget = ui->processingGlWidget;
    if (processingglWidget) {
        connect(this, &phaseEstwindow::setCustomScaleStatus, processingglWidget, &ProcessingGlWidget::setCustomScaleStatus);
        connect(this, &phaseEstwindow::setCustomScaleMin, processingglWidget, &ProcessingGlWidget::setCustomScaleMin);
        connect(this, &phaseEstwindow::setCustomScaleMax, processingglWidget, &ProcessingGlWidget::setCustomScaleMax);

        connect(this, &phaseEstwindow::getCustomScaleStatus, processingglWidget, &ProcessingGlWidget::getCustomScaleStatus);
        connect(this, &phaseEstwindow::getCustomScaleMin, processingglWidget, &ProcessingGlWidget::getCustomScaleMin);
        connect(this, &phaseEstwindow::getCustomScaleMax, processingglWidget, &ProcessingGlWidget::getCustomScaleMax);
        connect(this, &phaseEstwindow::setShowTriggers_A, processingglWidget, &ProcessingGlWidget::setShowTriggers_A);
        connect(this, &phaseEstwindow::setShowTriggers_B, processingglWidget, &ProcessingGlWidget::setShowTriggers_B);
        connect(this, &phaseEstwindow::setShowTriggers_out, processingglWidget, &ProcessingGlWidget::setShowTriggers_out);
        connect(this, &phaseEstwindow::switchPause, processingglWidget, &ProcessingGlWidget::switchPause);

        connect(this, &phaseEstwindow::updateWidgetChannelNames, processingglWidget, &ProcessingGlWidget::updateChannelNamesSTD);
        connect(this, &phaseEstwindow::setDrawXaxis, processingglWidget, &ProcessingGlWidget::setDrawXaxis);
        connect(this, &phaseEstwindow::updateTLineSpacing, processingglWidget, &ProcessingGlWidget::updateTLineSpacing);

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

    // Create the dock widget
    dockWidget = new QDockWidget("Phase error", this);
    dockWidget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    // Create the OpenGL widget
    HistogramWidget = new PolarHistogramOpenGLWidget(dockWidget);
    HistogramWidget->setFixedWidth(300);

    HistogramWidget->setMaximumHeight(700);

    // Set the OpenGL widget as the central widget of the dock
    dockWidget->setWidget(HistogramWidget);

    // Add the dock widget to the main window, but hide it by default
    addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    dockWidget->hide();
}

void phaseEstwindow::newEstStates(phaseEstimateStates states) {
    ui->checkBox_phaseEstimate->setChecked(states.performPhaseEstimation);
    ui->Filter_checkbox->setChecked(states.performFiltering);
    ui->checkBox_SNRcheck->setChecked(states.performSNRcheck);
    ui->checkBox_PhaseTargeting->setChecked(states.performPhaseTargeting);
    ui->checkBox_phaseError->setChecked(states.performPhaseDifference);
    ui->checkBox_Channels->setChecked(states.phasEst_display_all_EEG_channels);
}

void phaseEstwindow::updateSpatialChannelNames(std::vector<std::string> names) 
{
    spatial_channel_names = names;
}

phaseEstwindow::~phaseEstwindow()
{
    emit setphaseEstimateState(false);
    delete ui;
}

void phaseEstwindow::on_edge_editingFinished()
{
    bool ok;
    int value = ui->edge->text().toInt(&ok);
    if (ok) {
        phaseEstParams.edge = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void phaseEstwindow::on_modelOrder_editingFinished()
{
    bool ok;
    int value = ui->modelOrder->text().toInt(&ok);
    if (ok) {
        phaseEstParams.modelOrder = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void phaseEstwindow::on_hilbertLength_editingFinished()
{
    bool ok;
    int value = ui->hilbertLength->text().toInt(&ok);
    if (ok) {
        phaseEstParams.hilbertWinLength = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void phaseEstwindow::on_stimulationTarget_editingFinished()
{
    bool ok;
    double value = ui->stimulationTarget->text().toDouble(&ok);
    if (ok) {
        phaseEstParams.stimulation_target = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void phaseEstwindow::on_phaseShift_editingFinished()
{
    bool ok;
    int value = ui->phaseShift->text().toInt(&ok);
    if (ok) {
        phaseEstParams.phase_shift = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}




void phaseEstwindow::on_scaleMax_editingFinished()
{
    bool ok;
    double value = ui->scaleMax->text().toDouble(&ok);
    if (ok) {
        emit setCustomScaleMax(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void phaseEstwindow::on_scaleMin_editingFinished()
{
    bool ok;
    double value = ui->scaleMin->text().toDouble(&ok);
    if (ok) {
        emit setCustomScaleMin(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void phaseEstwindow::on_checkBox_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setCustomScaleStatus(isChecked);
}


void phaseEstwindow::on_Filter_checkbox_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setFilterState(isChecked);
}


void phaseEstwindow::on_checkBox_Channels_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setEEGViewState(isChecked);
}


void phaseEstwindow::on_checkBox_PhaseTargeting_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setPhaseTargetingState(isChecked);
}

void phaseEstwindow::on_checkBox_Stimulation_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    handler.setTriggerEnableStatus(isChecked);
}

void phaseEstwindow::on_checkBox_phaseEstimate_stateChanged(int arg1)
{
    if (handler.channelNamesSet()) {
        sampling_rate_ = handler.getSamplingRate();
        bool isChecked = (arg1 == Qt::Checked);

        //Begin Phase Estimation processing thread
        // emit startPhaseEstimationprocessing(phaseEstParams);
        emit setphaseEstimateState(isChecked);
        
        ui->comboBox_spatialTarget->clear(); // Clear the combo box before updating it

        for (const auto &name : spatial_channel_names)
        {
            ui->comboBox_spatialTarget->addItem(QString::fromStdString(name)); // Convert std::string to QString and add to combo box
        }

        setupComboBox_OuterElectrode(0);
    
    } else {
        ui->checkBox_phaseEstimate->setChecked(false);
        QMessageBox::warning(this, "Channel Naming Error", "Please load channel map in the EEG window.");
    }
}


void phaseEstwindow::on_setParamsButton_clicked()
{
    emit setPhaseEstParams(phaseEstParams);
}

void phaseEstwindow::on_comboBox_spatialTarget_currentIndexChanged(int index)
{
    emit setSpatilaTargetChannel(index);
    setupComboBox_OuterElectrode(index);
}

void phaseEstwindow::setupComboBox_OuterElectrode(int spatial_target_index) {

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
    connect(model, &QStandardItemModel::itemChanged, this, &phaseEstwindow::handleCheckboxChange);
}

void phaseEstwindow::handleCheckboxChange(QStandardItem* item) {
    int row = item->row();
    bool isChecked = item->checkState() == Qt::Checked;
    outerElectrodeCheckStates_[row] = isChecked;
    
    emit outerElectrodesStateChanged(outerElectrodeCheckStates_);
}

void phaseEstwindow::on_pushButton_pause_view_clicked()
{
    emit switchPause();
}

void phaseEstwindow::on_checkBox_triggers_A_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_A(isChecked);
}

void phaseEstwindow::on_checkBox_triggers_B_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_B(isChecked);
}

void phaseEstwindow::on_checkBox_phaseError_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setPhaseError(isChecked);
}

void phaseEstwindow::on_checkBox_Xaxis_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setDrawXaxis(isChecked);
}

void phaseEstwindow::on_lineEdit_XaxisSpacing_editingFinished()
{
    bool ok;
    int value = ui->lineEdit_XaxisSpacing->text().toInt(&ok);
    if (ok) {
        emit updateTLineSpacing(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}


void phaseEstwindow::on_checkBox_triggers_out_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setShowTriggers_out(isChecked);
}


void phaseEstwindow::on_comboBox_phaseError_currentIndexChanged(int index)
{
    emit setPhaseErrorType(index);
}


void phaseEstwindow::on_checkBox_SNRcheck_stateChanged(int arg1)
{
    bool isChecked = (arg1 == Qt::Checked);
    emit setSNRcheck(isChecked);
    HistogramWidget->clearSamples();
    showOpenGLDock();
}


void phaseEstwindow::on_SNRtreshold_editingFinished()
{
    bool ok;
    double value = ui->SNRtreshold->text().toDouble(&ok);
    if (ok) {
        emit setSNRthreshold(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}

void phaseEstwindow::newSNRmax_list(const std::vector<double>& list)
{
    // Initialize an empty QString to hold the tooltip text
    QString tooltipText;

    // Loop through the list and append each value to the tooltip text
    for (size_t i = 0; i < list.size(); ++i) {
        // Format each value with a label or index if desired
        tooltipText += QString("SNR max %1: %2").arg(i + 1).arg(list[i], 0, 'f', 2);

        // Add a newline character after each value except the last one
        if (i < list.size() - 1) {
            tooltipText += "\n";
        }
    }

    double Total_max = *std::max_element(list.begin(), list.end());

    tooltipText += "\n";
    tooltipText += QString("Max: %2").arg(Total_max, 0, 'f', 2);
    tooltipText += "\n";
    tooltipText += QString("Mean: %2").arg(std::accumulate(list.begin(), list.end(), 0.0) / list.size(), 0, 'f', 2);
    
    ui->SNRmax->setToolTip(tooltipText);
}


void phaseEstwindow::newSNRmax(double value) 
{
    ui->SNRmax->setText(QString::number(value));
}

void phaseEstwindow::on_SNRmax_editingFinished()
{
    bool ok;
    double value = ui->SNRmax->text().toDouble(&ok);
    if (ok) {
        emit sendSNRmax(value);
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number.");
    }
}

