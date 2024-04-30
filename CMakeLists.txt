cmake_minimum_required(VERSION 3.0)
project(MyMatplotlibCppProject)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Find PythonLibs (Preferably Python 3)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Include numpy
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
add_executable( real_time_eeg main.cpp 
                eeg_bridge/eeg_bridge.cpp 
                eeg_bridge/measurementStartPacket.cpp 
                eeg_bridge/samplePacket.cpp 
                eeg_bridge/networkUtils.cpp 
                dataHandler/dataHandler.cpp 
                dataHandler/GACorrection.cpp 
                dataProcessor/dataProcessor.cpp
                )

# Use the flag to suppress deprecated declarations warnings
target_compile_options(real_time_eeg PRIVATE -Wno-deprecated-declarations)

# Link Python libraries with keyword PRIVATE
target_link_libraries(real_time_eeg PRIVATE ${PYTHON_LIBRARIES})

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
