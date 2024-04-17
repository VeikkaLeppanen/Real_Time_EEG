#include "worker.h"
#include <iostream>

Worker::Worker(EegBridge &bridge, dataHandler &handler, volatile std::sig_atomic_t &signal_received, QObject* parent)
    : QObject(parent), bridge(bridge), handler(handler), signal_received(signal_received)
{ }

Worker::~Worker()
{ }

void Worker::process()
{
    // Ensure this code is thread-safe and does not interfere with the GUI thread
    try {
        bridge.bind_socket();
        bridge.spin(handler, signal_received);
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred: %1").arg(e.what()));
    }
}
