#include <vector>
#include <thread>
#include <mutex>
#include <omp.h>
#include <string>

#include "matplotlibcpp.h"
#include "dataHandler.h"

namespace plt = matplotlibcpp;
const uint8_t CHANNEL_COUNT = 3;
const uint32_t SAMPLING_RATE = 5000;
const uint32_t DOWNSAMPLING_FACTOR = 1;    // Down sampling currently not working
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
    // int channel_index = 0;
    int channel_count = CHANNEL_COUNT;

    plt::ion(); // Enable interactive mode
    plt::figure_size(1920, 100 * channel_count);

    // X values of the plot
    Eigen::VectorXd xVec = Eigen::VectorXd::LinSpaced(dataPoints, 0, DATABUFFER_LENGTH_IN_SECONDS);
    std::vector<double> x(xVec.data(), xVec.data() + xVec.size());

    // Y values of the plot
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(dataPoints);

    Eigen::MatrixXd bufferCopy(dataBuffer.rows(), dataPoints);
    while (true) { // Adjust this condition as needed
        plt::clf();

        {
            std::lock_guard<std::mutex> guard(dataMutex);
            bufferCopy = dataBuffer.getDataInOrder();
            // yVec = dataBuffer.getChannelDataInOrder(channel_index, DOWNSAMPLING_FACTOR);
        }
        
        for (int channel_index = 0; channel_index < channel_count; ++channel_index) {
            yVec = bufferCopy.row(channel_index);
            std::cout << yVec.size() << ' ' << x.size() << '\n';
            plt::subplot(channel_count, 1, channel_index + 1);
            plt::plot(x, std::vector<double>(yVec.data(), yVec.data() + yVec.size()));
            plt::title("Channel " + std::to_string(channel_index + 1)); // Optional: Add title to each subplot
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
