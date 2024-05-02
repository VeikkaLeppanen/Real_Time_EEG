#include "eegwindow.h"

#include <QLabel>

eegWindow::eegWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::EegWindow) {
    ui->setupUi(this);
    // QLabel *label = new QLabel("Hello from the second window!", this);
    setWindowTitle("Secondary Window");
    // resize(200, 100); // Set a default size
}

eegWindow::~eegWindow() {
    delete ui;  // Clean up the ui pointer
}
