cmake_minimum_required(VERSION 3.5)
project(Real_time_EEG VERSION 0.1 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# Find packages
find_package(QT NAMES Qt6 REQUIRED COMPONENTS Widgets OpenGL OpenGLWidgets) # Charts)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets OpenGL OpenGLWidgets) # Charts)
find_package(OpenGL REQUIRED)  # Added OpenGL package
find_package(Threads REQUIRED)
find_package(Boost 1.65 REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)

include_directories(/home/veikka/Work/Real_Time_EEG)
# include_directories(/usr/include/eigen3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/UI)
include_directories(${Boost_INCLUDE_DIRS})
# Add the include directories
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCES
    main.cpp
    eeg_bridge/eeg_bridge.cpp
    eeg_bridge/measurementStartPacket.cpp
    eeg_bridge/samplePacket.cpp
    eeg_bridge/networkUtils.cpp
    dataHandler/dataHandler.cpp
    dataHandler/GACorrection.cpp 
    dataProcessor/dataProcessor.cpp
    dataProcessor/processingFunctions.cpp
    UI/mainwindow.cpp
    UI/mainwindow.h
    UI/mainwindow.ui
    UI/eegwindow.cpp
    UI/eegwindow.h
    UI/eegwindow.ui
    UI/processingwindow.h 
    UI/processingwindow.cpp 
    UI/triggeringwindow.h 
    UI/triggeringwindow.cpp 
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

add_executable(real_time_eeg ${SOURCES}
    UI/triggeringwindow.h UI/triggeringwindow.cpp UI/triggeringwindow.ui)

target_link_libraries(real_time_eeg PRIVATE
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::OpenGL
    Qt${QT_VERSION_MAJOR}::OpenGLWidgets
    OpenGL::GL  # Linking OpenGL
    Threads::Threads
    ${Boost_LIBRARIES}
    fftw3
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(real_time_eeg)
endif()

# Use the flag to suppress deprecated declarations warnings
target_compile_options(real_time_eeg PRIVATE -Wno-deprecated-declarations) # -DEIGEN_USE_MKL_ALL

# Compiler optimization and OpenMP
include(CheckCXXCompilerFlag)
# Check and enable FMA instruction set if available
CHECK_CXX_COMPILER_FLAG("-mfma" COMPILER_SUPPORTS_MFMA)
if(COMPILER_SUPPORTS_MFMA)
    target_compile_options(real_time_eeg PRIVATE -mfma)
endif()

# Add compiler and linker options for OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_compile_options(real_time_eeg PRIVATE ${OpenMP_CXX_FLAGS})
    # Linking OpenMP libraries
    target_link_libraries(real_time_eeg PRIVATE ${OpenMP_CXX_LIBRARIES})
endif()

# Set other compiler optimizations
target_compile_options(real_time_eeg PRIVATE -O3)

