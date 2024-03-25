#include <vector>
#include <thread>
#include <mutex>
#include <omp.h>
#include <string>

#include "dataHandler/dataHandler.h"
#include "eeg_bridge/eeg_bridge.h"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

const uint8_t CHANNEL_COUNT = 4;
const uint32_t SAMPLING_RATE = 1000;
const uint32_t DELIVERY_RATE = 1000;
const uint32_t DOWNSAMPLING_FACTOR = 1;


// Data simulation loop
void dataSimulationLoop(dataHandler &handler) {
    
    std::cout << "Reseting dataHandler" << '\n';
    handler.reset_handler(CHANNEL_COUNT, SAMPLING_RATE, DELIVERY_RATE);
    
    std::cout << "Simulating data" << '\n';
    if (!handler.simulateData_sin()) {
        std::cout << "Simulation loop terminated" << '\n';
    }
}


// Loop for the data collecting and storing
void dataAcquisitionLoop(dataHandler &handler) {

    EegBridge bridge;
    bridge.bind_socket();
    bridge.spin(handler);
    
}

// Loop for matplotlibcpp graph visualization
void plottingLoop(dataHandler &handler) {
    
    // Wait for the handler initialization with measurementStartPackage
    while (!handler.isReady()) {
        // Optionally, sleep for a short duration to avoid busy waiting
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    int datapoints = handler.get_buffer_capacity();
    int dataPoints_downsampled = (datapoints + DOWNSAMPLING_FACTOR - 1) / DOWNSAMPLING_FACTOR;

    std::vector<int> channels_to_display = {0};
    int channel_count = CHANNEL_COUNT;

    plt::ion(); // Enable interactive mode
    plt::figure_size(1920, 100 * (channels_to_display.size() + 2));

    // X values of the plot
    Eigen::VectorXd xVec = Eigen::VectorXd::LinSpaced(dataPoints_downsampled, handler.get_buffer_length_in_seconds(), 0);
    std::vector<double> x(xVec.data(), xVec.data() + xVec.size());

    // Y values of the plot
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(dataPoints_downsampled);
    Eigen::MatrixXd downSampledData = Eigen::MatrixXd::Zero(handler.get_channel_count(), dataPoints_downsampled);

    while (true) { // Adjust this condition as needed

        plt::clf();

        int graph_ind = 0;
        for (int channel_index : channels_to_display) {

            yVec = handler.getChannelDataInOrder(channel_index, DOWNSAMPLING_FACTOR);

            // std::cout << "yVec: "; 
            // for (int abc = 0; abc < yVec.size(); abc++) {
            //     std::cout << yVec[abc] << ' ';
            // }
            // std::cout << '\n';

            // std::cout << yVec.size() << ' ' << x.size() << '\n';
            plt::subplot(channels_to_display.size() + 2, 1, graph_ind + 1);
            plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
            plt::title("Channel " + std::to_string(channel_index + 1)); // Optional: Add title to each subplot
            graph_ind++;
        }
        
        // Timestamp plotting
        // yVec = handler.getTimeStampsInOrder(DOWNSAMPLING_FACTOR);
        // plt::subplot(channels_to_display.size() + 2, 1, graph_ind + 1);
        // plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
        // plt::title("Time stamps");
        // graph_ind++;

        // Trigger plotting
        yVec = handler.getTriggersInOrder(DOWNSAMPLING_FACTOR);
        plt::subplot(channels_to_display.size() + 1, 1, graph_ind + 1);
        plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
        plt::title("Triggers");
        
        plt::pause(0.01); // Pause for a short period to allow the plot to update
    }
}

int main() {
    
    dataHandler handler;

    // std::thread dataThread(dataSimulationLoop, std::ref(handler));
    std::thread dataThread(dataAcquisitionLoop, std::ref(handler));
    std::this_thread::sleep_for(std::chrono::seconds(1));

    std::thread plotThread(plottingLoop, std::ref(handler));

    dataThread.join();
    plotThread.join();

    return 0;
}
