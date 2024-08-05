#include <vector>
#include <thread>
#include <mutex>
#include <string>
#include <csignal>
#include <chrono>
#include <fstream>
#include <malloc.h>
#include <sys/mman.h> // For mlockall

#include "dataProcessor/dataProcessor.h"
#include "dataProcessor/processingFunctions.h"
#include "dataHandler/dataHandler.h"
#include "eeg_bridge/eeg_bridge.h"

#include "UI/mainwindow.h"
#include <QApplication>
#include <QWidget>

// #include "matplotlibcpp.h"
// namespace plt = matplotlibcpp;

// In case of bind failed the previous process can be terminated with the following commands on linux
// lsof -i :50000       Find PID
// kill -9 PID          Kill the process

// This is used to terminate the program with Ctrl+C
volatile std::sig_atomic_t signal_received = 0;
void signal_handler(int signal) {
    signal_received = 1;
}

int main(int argc, char *argv[])
{
    // Lock all current and future memory pages into RAM
    if (mlockall(MCL_CURRENT | MCL_FUTURE) != 0) {
        std::cerr << "Failed to lock memory" << std::endl;
        return -1; // or handle the error appropriately
    }

    // Disable memory trimming
    if (mallopt(M_TRIM_THRESHOLD, -1) != 1) {
        std::cerr << "Failed to set M_TRIM_THRESHOLD" << std::endl;
    }

    // Prevent mmap from being used for memory allocation
    if (mallopt(M_MMAP_MAX, 0) != 1) {
        std::cerr << "Failed to set M_MMAP_MAX" << std::endl;
    }

    dataHandler handler;

    QApplication a(argc, argv);
    // qRegisterMetaType<Eigen::MatrixXd>("Eigen::MatrixXd");
    // qRegisterMetaType<Eigen::MatrixXd&>("Eigen::MatrixXd&");
    MainWindow w(handler, signal_received);
    w.show();
    return a.exec();
}

