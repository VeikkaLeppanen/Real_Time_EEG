#include "worker.h"

Worker::Worker(EegBridge &bridge, dataHandler &handler, volatile std::sig_atomic_t &signal_received, QObject* parent)
    : QObject(parent), bridge(bridge), handler(handler), signal_received(signal_received)
{ }

Worker::~Worker()
{ }

void Worker::process()
{
    // Set the current thread to use SCHED_RR
    pthread_t this_thread = pthread_self();
    struct sched_param params;
    params.sched_priority = sched_get_priority_max(SCHED_RR);
    
    if (pthread_setschedparam(this_thread, SCHED_RR, &params) != 0) {
        qDebug("Failed to set thread to real-time");
    } else {
        qDebug("Thread set to real-time successfully");
    }
    
    // Ensure this code is thread-safe and does not interfere with the GUI thread
    try {
        bridge.bind_socket();
        bridge.spin(handler, signal_received);
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred: %1").arg(e.what()));
    }
}
