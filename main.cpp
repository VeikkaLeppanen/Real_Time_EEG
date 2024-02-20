#include <vector>
#include <thread>
#include <mutex>
#include <omp.h>
#include <string>

#include "dataHandler.h"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

const uint8_t CHANNEL_COUNT = 20;
const uint32_t SAMPLING_RATE = 5000;
const uint32_t DELIVERY_RATE = 5000;
const uint32_t DOWNSAMPLING_FACTOR = 10;
const uint32_t SHORT_BUFFER_LENGTH_IN_SECONDS = 2;
const uint32_t LONG_BUFFER_LENGTH_IN_SECONDS = 10;

// Loop for the data collecting and storing
void dataAcquisitionLoop(circularEigenBuffer &short_data_buffer, circularEigenBuffer &long_data_buffer, std::mutex &dataMutex) {
    
    dataHandler handler(dataMutex, 
                        short_data_buffer, 
                        long_data_buffer, 
                        CHANNEL_COUNT,
                        SAMPLING_RATE,
                        DELIVERY_RATE);
    
    if (!handler.simulateData()) {
        std::cout << "Simulation loop terminated" << '\n';
    }
}

// Loop for matplotlibcpp graph visualization
void plottingLoop(circularEigenBuffer &short_data_buffer, std::mutex &dataMutex) {
    const int dataPoints = (SAMPLING_RATE * SHORT_BUFFER_LENGTH_IN_SECONDS + DOWNSAMPLING_FACTOR - 1) / DOWNSAMPLING_FACTOR;
    std::vector<int> channels_to_display = {0}; //, 5, 15};
    int channel_count = CHANNEL_COUNT;

    plt::ion(); // Enable interactive mode
    plt::figure_size(1920, 100 * channels_to_display.size());

    // X values of the plot
    Eigen::VectorXd xVec = Eigen::VectorXd::LinSpaced(dataPoints, 0, SHORT_BUFFER_LENGTH_IN_SECONDS);
    std::vector<double> x(xVec.data(), xVec.data() + xVec.size());

    // Y values of the plot
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(dataPoints);
    Eigen::MatrixXd downSampledData = Eigen::MatrixXd::Zero(short_data_buffer.getNumberOfRows(), dataPoints);

    while (true) { // Adjust this condition as needed
        plt::clf();

        int graph_ind = 0;
        for (int channel_index : channels_to_display) {
            {
            std::lock_guard<std::mutex> guard(dataMutex);
            yVec = short_data_buffer.getChannelDataInOrder(channel_index, DOWNSAMPLING_FACTOR);
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
    // Shared array of circularEigenBuffer to store data from channels. +1 for the time_stamps
    circularEigenBuffer short_data_buffer(CHANNEL_COUNT + 1, SAMPLING_RATE * SHORT_BUFFER_LENGTH_IN_SECONDS);
    circularEigenBuffer long_data_buffer(CHANNEL_COUNT + 1, SAMPLING_RATE * LONG_BUFFER_LENGTH_IN_SECONDS);
    // Mutex to protect the shared vector
    std::mutex dataMutex;

    std::thread dataThread(dataAcquisitionLoop, std::ref(short_data_buffer), std::ref(long_data_buffer), std::ref(dataMutex));
    std::thread plotThread(plottingLoop, std::ref(short_data_buffer), std::ref(dataMutex));

    dataThread.join();
    plotThread.join();

    return 0;
}
