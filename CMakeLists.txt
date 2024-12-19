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
find_package(OpenBLAS REQUIRED)

include_directories(/home/veikka/Work/Real_Time_EEG)
include_directories(/usr/local/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/UI)
include_directories(${Boost_INCLUDE_DIRS})
# Add the include directories
include_directories(${EIGEN3_INCLUDE_DIR} /usr/include/mkl)

if(OpenBLAS_FOUND)
    include_directories(${OpenBLAS_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "OpenBLAS not found")
endif()

file(GLOB MAIN_SOURCES "*.cpp" "*.h")
file(GLOB_RECURSE EEG_SOURCES "EEG/*.cpp" "EEG/*.h")
file(GLOB_RECURSE EEG_BRIDGE_SOURCES "devices/EEG/eeg_bridge/*.cpp" "devices/EEG/eeg_bridge/*.h")
file(GLOB_RECURSE DATAHANDLER_SOURCES "dataHandler/*.cpp" "dataHandler/*.h")
file(GLOB_RECURSE TMS_SOURCES "devices/TMS/*.cpp" "devices/TMS/*.h")
file(GLOB_RECURSE MATH_SOURCES "math/*.cpp" "math/*.h")
file(GLOB_RECURSE UTILS_SOURCES "utils/*.cpp" "utils/*.h")
file(GLOB_RECURSE UI_SOURCES "UI/*.cpp" "UI/*.h" "UI/*.ui")

# Combine all source lists
set(SOURCES ${MAIN_SOURCES} ${EEG_SOURCES} ${EEG_BRIDGE_SOURCES} ${DATAHANDLER_SOURCES} ${TMS_SOURCES} ${MATH_SOURCES} ${UTILS_SOURCES} ${UI_SOURCES})

add_definitions(-DEIGEN_USE_MKL_ALL)
add_definitions(-DEIGEN_USE_LAPACKE=0)

add_executable(real_time_eeg ${SOURCES})

# Link libraries
target_link_libraries(real_time_eeg PRIVATE
    Qt${QT_VERSION_MAJOR}::Widgets
    Qt${QT_VERSION_MAJOR}::OpenGL
    Qt${QT_VERSION_MAJOR}::OpenGLWidgets
    ${OpenBLAS_LIBRARIES}
    OpenGL::GL
    Threads::Threads
    ${Boost_LIBRARIES}
    fftw3
    LabJackM
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
