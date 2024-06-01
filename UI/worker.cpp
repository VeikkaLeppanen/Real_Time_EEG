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

    // Isolate the thread to core 0
    set_thread_affinity();
    
    // Ensure this code is thread-safe and does not interfere with the GUI thread
    try {
        bridge.bind_socket();
        bridge.spin(handler, signal_received);
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred: %1").arg(e.what()));
    }
}

void Worker::set_thread_affinity() {
    cpu_set_t cpuset;
    pthread_t thread = pthread_self();

    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);  // Assuming core 0 is reserved for the packet receiving thread

    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
}