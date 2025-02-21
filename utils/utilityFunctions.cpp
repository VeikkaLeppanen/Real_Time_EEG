#include "utilityFunctions.h"

// Function for writing Eigen matrix to CSV file
void writeMatrixdToCSV(const std::string& filename, const Eigen::MatrixXd& matrix) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j);
                if (j + 1 < matrix.cols()) file << ","; // Comma for next column
            }
            file << "\n"; // Newline for next row
        }
        file.close();
        std::cout << "Data written to " << filename << " successfully." << std::endl;
    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }
}

// Function for writing Eigen matrix to CSV file
void writeMatrixiToCSV(const std::string& filename, const Eigen::MatrixXi& matrix) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << matrix(i, j);
                if (j + 1 < matrix.cols()) file << ","; // Comma for next column
            }
            file << "\n"; // Newline for next row
        }
        file.close();
        std::cout << "Data written to " << filename << " successfully." << std::endl;
    } else {
        std::cerr << "Failed to open the file for writing." << std::endl;
    }
}

// Function for reading matrix form a CSV file
Eigen::MatrixXd readCSV(const std::string &file_path) {
    std::ifstream file(file_path);

    // Check if the file stream is valid
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::vector<double> matrix_entries;
    std::string line;
    int num_rows = 0;
    int num_cols = 0;

    // Read each line from the CSV file
    while (getline(file, line)) {
        std::stringstream line_stream(line);
        std::string cell;
        int current_row_cols = 0;

        // Split the line into comma-delimited elements
        while (getline(line_stream, cell, ',')) {
            matrix_entries.push_back(std::stod(cell)); // Convert string to double
            ++current_row_cols;
        }

        if (num_rows == 0) {
            num_cols = current_row_cols; // Set number of columns based on the first row
        }

        ++num_rows;
    }

    file.close();

    // Create an Eigen matrix from the vector of values
    Eigen::MatrixXd matrix(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            matrix(i, j) = matrix_entries[i * num_cols + j];
        }
    }

    return matrix;
}

// Function to convert Eigen::VectorXd to Eigen::MatrixXd for CSV writing
Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vec) {
    // Map the vector as a matrix with one column and vec.size() rows
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), vec.size(), 1);
}

// Function to convert std::vector<double> to a single-column Eigen::MatrixXd
Eigen::MatrixXd vectorToColumnMatrixd(const std::vector<double>& vec) {
    // Map the vector as a single-column matrix
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), vec.size(), 1);
}

// Function to convert std::vector<double> to a single-column Eigen::MatrixXd
Eigen::MatrixXi vectorToColumnMatrixi(const std::vector<int>& vec) {
    // Map the vector as a single-column matrix
    return Eigen::Map<const Eigen::MatrixXi>(vec.data(), vec.size(), 1);
}

// Function to convert std::vector<std::complex<double>> to Eigen::MatrixXd with real and imaginary parts in separate columns
Eigen::MatrixXd complexVectorToMatrix(const std::vector<std::complex<double>>& complexVec) {
    Eigen::MatrixXd combinedMatrix(complexVec.size(), 2);

    for (size_t i = 0; i < complexVec.size(); ++i) {
        combinedMatrix(i, 0) = complexVec[i].real();
        combinedMatrix(i, 1) = complexVec[i].imag();
    }

    return combinedMatrix;
}

// Function to calculate phase angle matrix from a complex vector
Eigen::MatrixXd phaseAngleToMatrix(const std::vector<std::complex<double>>& complexVec) {
    Eigen::MatrixXd phaseAngles(complexVec.size(), 1);

    for (size_t i = 0; i < complexVec.size(); ++i) {
        phaseAngles(i, 0) = std::arg(complexVec[i]);
    }

    return phaseAngles;
}