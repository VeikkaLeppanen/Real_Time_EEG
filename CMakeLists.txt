cmake_minimum_required(VERSION 3.5)
project(Real_time_EEG VERSION 0.1 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# Find packages
# find_package(PythonLibs 3 REQUIRED)
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

# Define QT_NO_KEYWORDS to avoid conflicts between Qt and Python
# add_definitions(-DQT_NO_KEYWORDS)

# Include directories
# include_directories(${PYTHON_INCLUDE_DIRS})
# execute_process(
#     COMMAND python3 -c "import numpy; print(numpy.get_include())"
#     OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(/home/veikka/Work/Real_Time_EEG)
include_directories(/usr/include/eigen3)

set(SOURCES
    main.cpp
    eeg_bridge/eeg_bridge.cpp
    eeg_bridge/measurementStartPacket.cpp
    eeg_bridge/samplePacket.cpp
    eeg_bridge/networkUtils.cpp
    dataHandler/dataHandler.cpp
    dataHandler/GACorrection.cpp
    UI/mainwindow.cpp
    UI/mainwindow.h
    UI/mainwindow.ui
    UI/worker.cpp
    UI/worker.h
)

add_executable(real_time_eeg ${SOURCES})

target_link_libraries(real_time_eeg PRIVATE Qt${QT_VERSION_MAJOR}::Widgets) # ${PYTHON_LIBRARIES}

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(real_time_eeg)
endif()

# Use the flag to suppress deprecated declarations warnings
target_compile_options(real_time_eeg PRIVATE -Wno-deprecated-declarations)
