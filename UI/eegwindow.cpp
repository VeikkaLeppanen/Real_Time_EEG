#include "eegwindow.h"

eegWindow::eegWindow(dataHandler &handler, volatile std::sig_atomic_t &signal_received, QWidget *parent) 
    : QMainWindow(parent), 
      ui(new Ui::EegWindow),
      handler(handler),
      signal_received(signal_received)
    {
        ui->setupUi(this);
        // QLabel *label = new QLabel("Hello from the second window!", this);
        setWindowTitle("EEG Window");
        // Optional: resize(800, 600); // Set a more appropriate default size

        Glwidget* glWidget = ui->openglWidget;  // Assuming the object name in the .ui is glWidget
        if (glWidget) {
            connect(glWidget, &Glwidget::fetchData, this, &eegWindow::updateData);
            connect(this, &eegWindow::updateChannelNames, glWidget, &Glwidget::updateChannelNames);
            
            emit updateChannelNames(handler.getChannelNames());
        } else {
            // Error handling if glWidget is not found
            qWarning("Glwidget not found in UI!");
        }
    }

eegWindow::~eegWindow() {
    delete ui;  // Clean up the ui pointer
}

void eegWindow::updateData() 
{
    Glwidget* glWidget = ui->openglWidget;  // Assuming ui->glWidget is the pointer to the Glwidget
    if (glWidget && handler.isReady()) {
        Eigen::MatrixXd newMatrix = handler.returnLatestDataInOrder(samples_to_display);
    
        glWidget->updateMatrix(newMatrix);  // You might need to implement updateMatrix in Glwidget
    }
}

void eegWindow::on_connectButton_clicked()
{
    emit connectEegBridge(port);

    while(!handler.isReady()) {
        QThread::msleep(500);
    }

    source_channels_ = handler.getSourceChannels();

    if (channelMap_.size() > 0) {

        std::vector<std::string> channelNames;

        for(size_t i = 0; i < source_channels_.size(); i++) {
            channelNames.push_back(channelMap_[source_channels_(i)]);
        }

        handler.setChannelNames(channelNames);

        emit updateChannelNames(channelNames);
    }
}

void eegWindow::handleError(const QString &error)
{
    // Error handling code here
    QMessageBox::critical(this, "Error", error);
}

void eegWindow::on_disconnectButton_clicked()
{
    signal_received = 1;
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

void eegWindow::on_GACorrectionStart_clicked()
{
    emit applyGACorrection(GALength, GAAverage);
}


void eegWindow::on_GACorrectionStop_clicked()
{
    emit stopGACorrection();
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

        std::vector<std::string> channelNames;

        for(size_t i = 0; i < source_channels_.size(); i++) {
            channelNames.push_back(channelMap_[source_channels_(i)]);
        }

        handler.setChannelNames(channelNames);

        emit updateChannelNames(channelNames);
    }
}

