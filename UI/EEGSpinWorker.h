#ifndef WORKER_H
#define WORKER_H

#include <QObject>
#include <iostream>
#include <pthread.h>
#include <sched.h>

#include "../dataHandler/dataHandler.h"
#include "../devices/EEG/eeg_bridge/eeg_bridge.h"

class EegBridge;
class dataHandler;

class EEGSpinWorker : public QObject {
    Q_OBJECT

public:
    explicit EEGSpinWorker(EegBridge &bridge, dataHandler &handler, volatile std::sig_atomic_t &signal_received, QObject* parent = nullptr);
    ~EEGSpinWorker();

signals:
    void finished();
    void error(QString err);

public slots:
    void process();

private:
    EegBridge &bridge;
    dataHandler &handler;
    volatile std::sig_atomic_t &signal_received;

    void bridge_handler_spin(EegBridge &bridge, dataHandler &handler, volatile std::sig_atomic_t &signal_received);

    void set_thread_affinity();
};

#endif // WORKER_H
