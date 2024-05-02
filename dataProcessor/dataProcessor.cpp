#include "dataProcessor.h"
// Downsampling function
Eigen::MatrixXd dataProcessor::downsample(const Eigen::MatrixXd& input, int factor) {
    if (factor <= 0) {
        throw std::invalid_argument("Downsampling factor must be greater than zero.");
    }

    // Compute the number of columns in the downsampled matrix
    int newCols = (input.cols() + factor - 1) / factor;

    // Initialize the downsampled matrix
    Eigen::MatrixXd output(input.rows(), newCols);

    for (int j = 0, col = 0; j < input.cols() && col < newCols; j += factor, ++col) {
        output.col(col) = input.col(j);
    }

    return output;
}

Eigen::MatrixXd dataProcessor::delayEmbed(const Eigen::MatrixXd& CWL, int delay) {
    int rows = CWL.rows();
    int cols = CWL.cols();
    Eigen::MatrixXd expandedCWL(rows, cols + 2 * delay);

    // Pad columns with zeros based on delay
    expandedCWL.block(0, delay, rows, cols) = CWL;

    // Fill in zeros for padding
    expandedCWL.block(0, 0, rows, delay).setZero();
    expandedCWL.block(0, cols + delay, rows, delay).setZero();

    return expandedCWL;
}

Eigen::MatrixXd dataProcessor::removeBCG(const Eigen::MatrixXd& EEG, Eigen::MatrixXd CWL, int delay) {
    if (delay > 0) {
        // Flip and perform delay embedding if delay is positive
        CWL = delayEmbed(CWL.colwise().reverse(), 1 + 2 * delay).colwise().reverse();
    }

    // Calculate pseudo-inverse of the expanded CWL matrix
    Eigen::MatrixXd inv_expCWL = CWL.completeOrthogonalDecomposition().pseudoInverse();

    // Square matrix for fits calculation
    Eigen::MatrixXd squareCWL = inv_expCWL * CWL;

    // Initialize fits matrix
    Eigen::MatrixXd fits = EEG * squareCWL;

    // Subtract fits from EEG to get BCG corrected data
    Eigen::MatrixXd data_BCG = EEG - fits;

    return data_BCG;
}

void dataProcessor::process(Eigen::MatrixXd EEG, Eigen::MatrixXd CWL) {
    std::cout << "Processing data on thread " << std::this_thread::get_id() << std::endl;
    
    // LOWPASS filter 120Hz   BANDPASS 0.33-125Hz


    // DONWSAMPLING
    int downsampling_factor = 10;
    Eigen::MatrixXd EEG_downSampled = downsample(EEG, downsampling_factor);
    Eigen::MatrixXd CWL_downSampled = downsample(CWL, downsampling_factor);


    // removeBCG
    int delay = 0;
    Eigen::MatrixXd correctedEEG = removeBCG(EEG, CWL, delay);

    std::cout << "Data process completed on thread " << std::this_thread::get_id() << std::endl;
}










// More efficient alternative to getMultipleChannelDataInOrder that can be used if the wanted channels are adjacent. Retrieves adjacent channels as a block
Eigen::MatrixXd dataProcessor::getBlockChannelDataInOrder(int first_channel_index, int number_of_channels, int number_of_samples) {

    Eigen::MatrixXd outputData(number_of_channels, number_of_samples);
    size_t channel_index = 0;


    // Calculate the number of samples that fit before reaching the end
    int fitToEnd = std::min(number_of_samples, (buffer_capacity_ - static_cast<int>(current_data_index_)));
    int overflow = number_of_samples - fitToEnd;
    
    if (fitToEnd > 0) {
        outputData.leftCols(fitToEnd) = sample_buffer_.block(first_channel_index, current_data_index_, number_of_channels, fitToEnd);
    }
    
    if (overflow > 0) {
        outputData.rightCols(overflow) = sample_buffer_.block(first_channel_index, 0, number_of_channels, overflow);
    }

    return outputData;
}