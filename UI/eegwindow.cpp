#include "eegwindow.h"

eegWindow::eegWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent) 
    : QMainWindow(parent), 
      ui(new Ui::EegWindow),
      handler(handler),
      signal_received(signal_received)
    {
        ui->setupUi(this);
        ui->lineEditPort->setValidator(new QIntValidator(0, 65535, this));
        ui->lineEditTimeOut->setValidator(new QIntValidator(0, 1000, this));
        ui->lineEditGALength->setValidator(new QIntValidator(0, 500000, this));
        ui->lineEditGAaverage->setValidator(new QIntValidator(0, 1000, this));
        setWindowTitle("EEG Window");
        // Optional: resize(800, 600); // Set a more appropriate default size

        Glwidget* glWidget = ui->openglWidget;
        if (glWidget) {
            connect(glWidget, &Glwidget::fetchData, this, &eegWindow::updateData);
            connect(this, &eegWindow::updateChannelNamesSTD, glWidget, &Glwidget::updateChannelNamesSTD);
            connect(this, &eegWindow::updateChannelNamesQt, glWidget, &Glwidget::updateChannelNamesQt);
            connect(this, &eegWindow::updateChannelDisplayState, glWidget, &Glwidget::updateChannelDisplayState);
            
            emit updateChannelNamesSTD(handler.getChannelNames());
        } else {
            // Error handling if glWidget is not found
            qWarning("Glwidget not found in UI!");
        }
    }

eegWindow::~eegWindow() {
    delete ui;
}

void eegWindow::handleError(const QString &error)
{
    // Error handling code here
    QMessageBox::critical(this, "Error", error);
}

void eegWindow::updateData() 
{
    Glwidget* glWidget = ui->openglWidget;
    if (glWidget && handler.isReady()) {
        Eigen::MatrixXd newMatrix = handler.returnLatestDataInOrder(samples_to_display);
    
        glWidget->updateMatrix(newMatrix);
    }
}

void eegWindow::on_connectButton_clicked()
{
    emit connectEegBridge(port, bridge_timeout);

    while(!handler.isReady()) {
        QThread::msleep(500);
    }

    source_channels_ = handler.getSourceChannels();

    if (channelMap_.size() > 0) {
        setupChannelNames();
    } else {
        QStringList QchannelNames;
        for(size_t i = 0; i < source_channels_.size(); i++) {
            QchannelNames.append("Undefined");
        }

        channelNames_ = QchannelNames;

        emit updateChannelNamesQt(QchannelNames);
        setupComboBox();
    }
}

void eegWindow::on_disconnectButton_clicked()
{
    signal_received = 1;
}

void eegWindow::on_sourceChannelLoad_clicked()
{
    std::vector<std::string> channelMapStd;
    QString fileName = QFileDialog::getOpenFileName(
        this,                 // parent widget
        "Open Document",      // dialog caption
        QDir::homePath(),     // starting directory
        "Text Files (*.txt)"  // file types
    );

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::ReadOnly)) {
            QTextStream in(&file);
            QString line = in.readLine();  // Read and discard the header line

            while (!in.atEnd()) {
                line = in.readLine();
                QStringList fields = line.split('\t');  // assuming tab-separated values
                if (fields.size() > 1) {  // Check if there is at least one field for input number and one for name
                    channelMapStd.push_back(fields[1].toStdString());
                }
            }
            file.close();
        } else {
            // Handle error if the file can't be opened
            qDebug("Failed to open the file for reading.");
        }
    }

    channelMap_ = channelMapStd;

    if (source_channels_.size() > 0) {
        setupChannelNames();
    }
}

void eegWindow::setupChannelNames()
{
    std::vector<std::string> channelNames;
    QStringList QchannelNames;
    for(size_t i = 0; i < source_channels_.size(); i++) {
        std::string C_Name = channelMap_[source_channels_(i) - 1];              // CHECK LATER FOR CORRECT INDEXING!
        channelNames.push_back(C_Name);
        QchannelNames.append(QString::fromStdString(C_Name));
    }
    handler.setChannelNames(channelNames);
    
    channelNames_ = QchannelNames;

    emit updateChannelNamesQt(QchannelNames);
    setupComboBox();
}

void eegWindow::setupComboBox() {
    QStandardItemModel* model = new QStandardItemModel(channelNames_.size(), 1, this);

    channelCheckStates_.resize(channelNames_.size(), true);

    for (int i = 0; i < channelNames_.size(); ++i) {
        QStandardItem* item = new QStandardItem(channelNames_.at(i));
        item->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
        item->setData(Qt::Checked, Qt::CheckStateRole);
        model->setItem(i, 0, item);
    }

    ui->comboBox->setModel(model);
    connect(model, &QStandardItemModel::itemChanged, this, &eegWindow::handleCheckboxChange);
    emit updateChannelDisplayState(channelCheckStates_);
}


void eegWindow::handleCheckboxChange(QStandardItem* item) {
    int row = item->row();
    bool isChecked = item->checkState() == Qt::Checked;
    channelCheckStates_[row] = isChecked;  // Directly set the boolean state

    // You might want to do something immediately after the change or just store the state.
    emit updateChannelDisplayState(channelCheckStates_);
}


void eegWindow::on_lineEditPort_editingFinished()
{
    bool ok;
    int value = ui->lineEditPort->text().toInt(&ok);
    if (ok) {
        port = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 65535.");
    }
}

void eegWindow::on_lineEditTimeOut_editingFinished()
{
    bool ok;
    int value = ui->lineEditTimeOut->text().toInt(&ok);
    if (ok) {
        bridge_timeout = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 1000.");
    }
}

void eegWindow::on_lineEditGraphSamples_editingFinished()
{
    bool ok;
    int value = ui->lineEditGraphSamples->text().toInt(&ok);
    int handler_capasity = handler.get_buffer_capacity();
    if (ok && (0 < value) && (value <= handler_capasity)) {
        samples_to_display = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 1 and " + QString::number(handler_capasity));
    }
}

void eegWindow::on_lineEditGALength_editingFinished()
{
    bool ok;
    int value = ui->lineEditGALength->text().toInt(&ok);
    if (ok) {
        GALength = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 500000.");
    }
}

void eegWindow::on_lineEditGAaverage_editingFinished()
{
    bool ok;
    int value = ui->lineEditGAaverage->text().toInt(&ok);
    if (ok) {
        GAAverage = value;
    } else {
        QMessageBox::warning(this, "Input Error", "Please enter a valid number between 0 and 1000.");
    }
}

void eegWindow::on_HandlerApplyButton_clicked()
{
    emit applyGACorrection(GALength, GAAverage);
}

void eegWindow::on_GACorrectionStart_clicked()
{
    emit startGACorrection();
}


void eegWindow::on_GACorrectionStop_clicked()
{
    emit stopGACorrection();
}