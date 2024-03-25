cmake_minimum_required(VERSION 3.0)
project(MyMatplotlibCppProject)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Find PythonLibs (Preferably Python 3)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Include numpy if you're using it
execute_process(
    COMMAND python3 -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
include_directories(${NUMPY_INCLUDE_DIRS})

# Include directory for matplotlibcpp and Eigen
include_directories(/home/veikka/Work/Real_Time_EEG)
include_directories(/usr/include/eigen3)

# Add executable
add_executable(my_project main.cpp eeg_bridge/eeg_bridge.cpp eeg_bridge/measurementStartPacket.cpp eeg_bridge/samplePacket.cpp eeg_bridge/networkUtils.cpp dataHandler/dataHandler.cpp dataHandler/GACorrection.cpp)

# Use the flag to suppress deprecated declarations warnings
target_compile_options(my_project PRIVATE -Wno-deprecated-declarations)

# Link Python libraries
target_link_libraries(my_project ${PYTHON_LIBRARIES})
