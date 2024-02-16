#include <vector>
#include <thread>
#include <mutex>
#include <omp.h>
#include <string>

#include "matplotlibcpp.h"
#include "dataHandler.h"

namespace plt = matplotlibcpp;
const uint8_t CHANNEL_COUNT = 20;
const uint32_t SAMPLING_RATE = 5000;
const uint32_t DOWNSAMPLING_FACTOR = 10;
const uint32_t DATABUFFER_LENGTH_IN_SECONDS = 5;

// Loop for the data collecting and storing
void dataAcquisitionLoop(circularEigenBuffer &dataBuffer, std::mutex &dataMutex) {
    
    dataHandler handler(dataBuffer, dataMutex, SAMPLING_RATE * DATABUFFER_LENGTH_IN_SECONDS, CHANNEL_COUNT);
    
    if (!handler.simulateData(SAMPLING_RATE)) {
        std::cout << "Simulation loop terminated" << '\n';
    }
}

void plottingLoop(circularEigenBuffer &dataBuffer, std::mutex &dataMutex) {
    const int dataPoints = SAMPLING_RATE * DATABUFFER_LENGTH_IN_SECONDS / DOWNSAMPLING_FACTOR;
    std::vector<int> channels_to_display = {1, 5, 15};
    int channel_count = CHANNEL_COUNT;

    plt::ion(); // Enable interactive mode
    plt::figure_size(1920, 100 * channels_to_display.size());

    // X values of the plot
    Eigen::VectorXd xVec = Eigen::VectorXd::LinSpaced(dataPoints, 0, DATABUFFER_LENGTH_IN_SECONDS);
    std::vector<double> x(xVec.data(), xVec.data() + xVec.size());

    // Y values of the plot
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(dataPoints);

    while (true) { // Adjust this condition as needed
        plt::clf();
        
        int graph_ind = 0;
        for (int channel_index : channels_to_display) {
            {
            std::lock_guard<std::mutex> guard(dataMutex);
            yVec = dataBuffer.getChannelDataInOrder(channel_index, DOWNSAMPLING_FACTOR);
            }

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
    // Shared array of circularEigenBuffer to store data from channels
    circularEigenBuffer dataBuffer(CHANNEL_COUNT, SAMPLING_RATE * DATABUFFER_LENGTH_IN_SECONDS);
    // Mutex to protect the shared vector
    std::mutex dataMutex;

    std::thread dataThread(dataAcquisitionLoop, std::ref(dataBuffer), std::ref(dataMutex));
    std::thread plotThread(plottingLoop, std::ref(dataBuffer), std::ref(dataMutex));

    dataThread.join();
    plotThread.join();

    return 0;
}
