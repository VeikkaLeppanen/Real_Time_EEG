#ifndef WORKER_H
#define WORKER_H

#include <QObject>
#include "../dataHandler/dataHandler.h"
#include "../eeg_bridge/eeg_bridge.h"

class EegBridge;
class dataHandler;

class Worker : public QObject {
    Q_OBJECT

public:
    explicit Worker(EegBridge &bridge, dataHandler &handler, volatile std::sig_atomic_t &signal_received, QObject* parent = nullptr);
    ~Worker();

signals:
    void finished();
    void error(QString err);

public slots:
    void process();

private:
    EegBridge &bridge;
    dataHandler &handler;
    volatile std::sig_atomic_t &signal_received;
};

#endif // WORKER_H
