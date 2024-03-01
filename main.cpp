#include <vector>
#include <thread>
#include <mutex>
#include <omp.h>
#include <string>

#include "dataHandler/dataHandler.h"
#include "eeg_bridge/eeg_bridge.h"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

const uint8_t CHANNEL_COUNT = 20;
const uint32_t SAMPLING_RATE = 5000;
const uint32_t DELIVERY_RATE = 5000;
const uint32_t DOWNSAMPLING_FACTOR = 10;
const uint32_t SHORT_BUFFER_LENGTH_IN_SECONDS = 2;
const uint32_t LONG_BUFFER_LENGTH_IN_SECONDS = 6;


// Data simulation loop
void dataSimulationLoop(dataHandler &handler) {
    
    handler.reset_handler(CHANNEL_COUNT, SAMPLING_RATE, DELIVERY_RATE);
    
    if (!handler.simulateData()) {
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
    int datapoints = handler.get_short_buffer_capacity();
    int dataPoints_downsampled = (datapoints + DOWNSAMPLING_FACTOR - 1) / DOWNSAMPLING_FACTOR;

    std::vector<int> channels_to_display = {0}; //, 20};
    int channel_count = CHANNEL_COUNT;

    plt::ion(); // Enable interactive mode
    plt::figure_size(1920, 100 * channels_to_display.size());

    // X values of the plot
    Eigen::VectorXd xVec = Eigen::VectorXd::LinSpaced(dataPoints_downsampled, 0, SHORT_BUFFER_LENGTH_IN_SECONDS);
    std::vector<double> x(xVec.data(), xVec.data() + xVec.size());

    // Y values of the plot
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(dataPoints_downsampled);
    Eigen::MatrixXd downSampledData = Eigen::MatrixXd::Zero(handler.get_short_buffer_num_rows(), dataPoints_downsampled);

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
            plt::subplot(channels_to_display.size(), 1, graph_ind + 1);
            plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
            plt::title("Channel " + std::to_string(channel_index + 1)); // Optional: Add title to each subplot
            graph_ind++;
        }
        
        plt::pause(0.01); // Pause for a short period to allow the plot to update
    }
}

int main() {
    
    dataHandler handler;

    std::thread dataThread(dataSimulationLoop, std::ref(handler));
    // std::thread dataThread(dataAcquisitionLoop, std::ref(handler));
    std::thread plotThread(plottingLoop, std::ref(handler));

    dataThread.join();
    plotThread.join();

    return 0;
}

/*
int main() {
    
    // // Shared circularEigenBuffers to store data from channels. +1 for the time_stamps
    circularEigenBuffer short_data_buffer(CHANNEL_COUNT + 1, SAMPLING_RATE * SHORT_BUFFER_LENGTH_IN_SECONDS);
    circularEigenBuffer long_data_buffer(CHANNEL_COUNT + 1, SAMPLING_RATE * LONG_BUFFER_LENGTH_IN_SECONDS);

    // // Mutex to protect the shared vector
    std::mutex dataMutex;

    std::thread dataThread(dataAcquisitionLoop_simulate, std::ref(short_data_buffer), std::ref(long_data_buffer), std::ref(dataMutex));
    std::thread plotThread(plottingLoop, std::ref(short_data_buffer), std::ref(dataMutex));

    dataThread.join();
    plotThread.join();

    return 0;
}
*/