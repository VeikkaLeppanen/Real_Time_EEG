cmake_minimum_required(VERSION 3.0)
project(MyMatplotlibCppProject)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Find PythonLibs and Qt6
find_package(PythonLibs 3 REQUIRED)
find_package(Qt6 COMPONENTS Widgets REQUIRED)

# Define QT_NO_KEYWORDS to avoid conflicts between Qt and Python
add_definitions(-DQT_NO_KEYWORDS)

# Include directories
include_directories(${PYTHON_INCLUDE_DIRS})
execute_process(
    COMMAND python3 -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(/home/veikka/Work/Real_Time_EEG)
include_directories(/usr/include/eigen3)

# Add executable
add_executable(real_time_eeg main.cpp eeg_bridge/eeg_bridge.cpp eeg_bridge/measurementStartPacket.cpp eeg_bridge/samplePacket.cpp eeg_bridge/networkUtils.cpp dataHandler/dataHandler.cpp dataHandler/GACorrection.cpp)

# Use the flag to suppress deprecated declarations warnings
target_compile_options(real_time_eeg PRIVATE -Wno-deprecated-declarations)

# Link Python and Qt6 libraries
target_link_libraries(real_time_eeg ${PYTHON_LIBRARIES} Qt6::Widgets)

# Enable Qt6 features
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
