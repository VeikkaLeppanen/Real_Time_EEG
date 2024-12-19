#include "EEGSpinWorker.h"

EEGSpinWorker::EEGSpinWorker(EegBridge &bridge, dataHandler &handler, volatile std::sig_atomic_t &signal_received, QObject* parent)
    : QObject(parent), bridge(bridge), handler(handler), signal_received(signal_received)
{ }

EEGSpinWorker::~EEGSpinWorker()
{ }

void EEGSpinWorker::process()
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
        // bridge.spin(handler, signal_received);
        bridge_handler_spin(bridge, handler, signal_received);
        emit finished();
    } catch (std::exception& e) {
        emit error(QString("An error occurred: %1").arg(e.what()));
    }
}

void EEGSpinWorker::set_thread_affinity() {
    cpu_set_t cpuset;
    pthread_t thread = pthread_self();

    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);  // Assuming core 0 is reserved for the packet receiving thread

    pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
}

void EEGSpinWorker::bridge_handler_spin(EegBridge &bridge, dataHandler &handler, volatile std::sig_atomic_t &signal_received) {
    
    try {
    bridge.running = true;
    std::cout << "Waiting for measurement start..." << '\n';
    bridge.eeg_bridge_status = WAITING_MEASUREMENT_START;
    while (!signal_received) {
        int n = bridge.receive_packet();
        if (n <= 0) {
            if (signal_received) break; // Check if the signal caused recvfrom to fail
            std::cerr << "Receive failed" << '\n';
            // break; // Optionally break on other errors too
            continue;
        }

        unsigned char firstByte = bridge.buffer[0];
        // Handle packets
        switch (firstByte)
        {
        case 0x01: { // MeasurementStartPacket
            // Add buffer validation
            if (n <= 0 || bridge.buffer == nullptr) {
                throw std::runtime_error("Invalid buffer or size in MeasurementStart: n=" + std::to_string(n));
            }
            
            std::cout << "MeasurementStart package received!" << '\n';
            std::cout << "Packet size: " << n << " Bytes" << '\n';
            measurement_start_packet packet_info;
            std::vector<uint16_t> SourceChannels;
            std::vector<uint8_t> ChannelTypes;
            
            try {
                deserializeMeasurementStartPacket_pointer(bridge.buffer, n, packet_info, SourceChannels, ChannelTypes);
            } catch (const std::exception& e) {
                std::cerr << "MeasurementStart deserialization error: " << e.what() << '\n';
                break;
            }

            // Validate channel data
            if (SourceChannels.empty()) {
                throw std::runtime_error("No channels received in MeasurementStart packet");
            }
            
            // Divide channels into data and trigger sources
            std::vector<uint16_t> data_channel_sources;
            uint16_t trigger_channel_source = -1;
            for (size_t i = 0; i < SourceChannels.size(); i++) {
                uint16_t source = SourceChannels[i];
                if(source < 60000) { 
                    data_channel_sources.push_back(source);
                } else {
                    trigger_channel_source = source; 
                }
            }

            bridge.numChannels = packet_info.NumChannels;
            
            if (trigger_channel_source != -1) bridge.numDataChannels = bridge.numChannels - 1;              // Excluding trigger channel
            else bridge.numDataChannels = bridge.numChannels;
            
            bridge.sampling_rate = packet_info.SamplingRateHz;
            bridge.lastSequenceNumber = -1;

            handler.setSourceChannels(data_channel_sources);
            handler.setTriggerSource(trigger_channel_source);

            bridge.data_handler_samples = Eigen::MatrixXd::Zero(bridge.numDataChannels, 10);

            std::cout << "MeasurementStart package processed!\n";

            handler.reset_handler(bridge.numDataChannels, bridge.sampling_rate);
            std::cout << "DataHandler reset!\n";
            bridge.eeg_bridge_status = MEASUREMENT_IN_PROGRESS;
            std::cout << "Waiting for packets..." << '\n';

            break;

        } case 0x02: { // SamplesPacket
            if (bridge.eeg_bridge_status == WAITING_MEASUREMENT_START || !handler.isReady()) break;

            // Add buffer validation
            if (n <= 0 || bridge.buffer == nullptr) {
                throw std::runtime_error("Invalid buffer or size: n=" + std::to_string(n));
            }

            // Deserialize the received data into a sample_packet instance
            sample_packet packet_info;
            Eigen::VectorXi triggers_A;
            Eigen::VectorXi triggers_B;

            try {
                deserializeSamplePacketEigen_pointer(bridge.buffer, n, packet_info, 
                    bridge.data_handler_samples, triggers_A, triggers_B, 
                    (bridge.numChannels > bridge.numDataChannels));
            } catch (const std::exception& e) {
                std::cerr << "Deserialization error: " << e.what() << '\n';
                break;
            }
            
            int sequenceNumber = packet_info.PacketSeqNo;

            if (bridge.numDataChannels != bridge.data_handler_samples.rows()) {
                std::cerr << "Error: numDataChannels is not equal to total rows in data_handler_samples. " << bridge.numDataChannels << ' ' << bridge.data_handler_samples.rows() << '\n';
                break;
            }

            Eigen::MatrixXd data_samples = ((bridge.data_handler_samples * bridge.DC_MODE_SCALE) / bridge.NANO_TO_MICRO_CONVERSION);

            // Ensure column access is valid
            if (packet_info.NumSampleBundles > data_samples.cols()) {
                std::cerr << "Error: Requested more sample bundles than available columns in data_samples." << '\n';
                break;
            }

            for (int i = 0; i < packet_info.NumSampleBundles; i++) {
                if (i >= data_samples.cols()) {
                    std::cerr << "Error: Column index " << i << " is out of range for data_samples with columns " << data_samples.cols() << '\n';
                    break;
                }
                handler.addData(data_samples.col(i), static_cast<double>(packet_info.FirstSampleTime), triggers_A(i), triggers_B(i), sequenceNumber);
            }

            // Check for dropped packets
            if (bridge.lastSequenceNumber != -1 && sequenceNumber != (bridge.lastSequenceNumber + 1)) {
                std::cerr << "Packet loss detected. Expected sequence: " << (bridge.lastSequenceNumber + 1) << ", but received: " << sequenceNumber << '\n';
                // TODO: Handle packet loss
            }

            bridge.lastSequenceNumber = sequenceNumber; // Update the latest sequence number

            // Debug output to confirm data integrity
            // std::cout << "Package " << sequenceNumber << " received!\n";
            // std::cout << "Channels: " << data_samples.rows() << ", " << data_samples.cols() << '\n';
            break;

        } case 0x03: { // TriggerPacket
            /* code */
            break;

        } case 0x04: { // MeasurementEndPacket
            /* code */
            break;
        
        } case 0x05: { // HardwareStatePacket
            /* code */
            break;
        
        
        } default:
            break;
        }
    }

    bridge.running = false;
    std::cout << "Shutting down..." << '\n';
    bridge.close_socket();

    } catch (const std::exception& e) {
        std::cerr << "Eeg_bridge exception: " << e.what() << '\n';
        std::cerr << boost::stacktrace::stacktrace();
    }
}