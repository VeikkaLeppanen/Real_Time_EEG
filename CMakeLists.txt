cmake_minimum_required(VERSION 3.5)
project(Real_time_EEG VERSION 0.1 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# Find packages
find_package(QT NAMES Qt6 REQUIRED COMPONENTS Widgets Charts)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets Charts)
find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.65 REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(/home/veikka/Work/Real_Time_EEG)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/UI)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${MKLROOT}/include)
include_directories(/usr/include/mkl)

# Sources
set(SOURCES
    main.cpp
    eeg_bridge/eeg_bridge.cpp
    eeg_bridge/measurementStartPacket.cpp
    eeg_bridge/samplePacket.cpp
    eeg_bridge/networkUtils.cpp
    dataHandler/dataHandler.cpp
    dataHandler/GACorrection.cpp 
# dataProcessor/dataProcessor.cpp
    dataProcessor/processingFunctions.cpp
    UI/mainwindow.cpp
    UI/mainwindow.h
    UI/mainwindow.ui
    UI/eegwindow.cpp
    UI/eegwindow.h
    UI/eegwindow.ui
    UI/processingwindow.h 
    UI/processingwindow.cpp 
    UI/processingwindow.ui
    UI/processingglwidget.cpp 
    UI/processingglwidget.h
    UI/processingworker.cpp
    UI/processingworker.h
    UI/worker.cpp
    UI/worker.h
    UI/glwidget.cpp
    UI/glwidget.h
    UI/mainglwidget.cpp
    UI/mainglwidget.h
)

add_definitions(-DEIGEN_USE_MKL_ALL)
add_definitions(-DEIGEN_USE_LAPACKE=0)

add_executable(real_time_eeg ${SOURCES})

# Link libraries
target_link_libraries(real_time_eeg PRIVATE
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::Charts
    OpenGL::GL
    Threads::Threads
    ${Boost_LIBRARIES}
    fftw3
    -m64 -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(real_time_eeg)
endif()

# Compiler flags for MKL and Eigen
target_compile_definitions(real_time_eeg PRIVATE)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-mfma" COMPILER_SUPPORTS_MFMA)
if(COMPILER_SUPPORTS_MFMA)
    target_compile_options(real_time_eeg PRIVATE -mfma)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_compile_options(real_time_eeg PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(real_time_eeg PRIVATE ${OpenMP_CXX_LIBRARIES})
endif()

target_compile_options(real_time_eeg PRIVATE -O3 -DMKL_ILP64 -m64 -I"${MKLROOT}/include")