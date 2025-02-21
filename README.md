# Real_Time_EEG

A real-time EEG data acquisition, processing, and visualization system with integrated TMS (Transcranial Magnetic Stimulation) control and MRI visualization capabilities.

## Features

- Real-time EEG data acquisition and processing
- TMS device control (MagPro) via serial communication
- MRI/fMRI data visualization and ROI analysis
- Real-time signal filtering and phase estimation
- Multi-threaded data handling with Qt-based GUI
- CSV data import/export capabilities

## Prerequisites

- Qt 6.7.0 or higher
- CMake 3.5 or higher
- Boost 1.65 or higher
- Eigen3
- FFTW3
- LabJackM
- Nibrary v0.3.0
- OpenGL
- C++17 compatible compiler

## Building

1. Set up Nibrary:
```bash
export NIBRARY_ROOT="/path/to/nibrary/build-static"
```

2. Create build directory and compile:
```bash
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=~/Qt/6.7.0/gcc_64/lib/cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

3. Set required permissions:
```bash
sudo setcap cap_ipc_lock+ep real_time_eeg
sudo chmod 666 /dev/ttyUSB0
```

## Usage

Run the compiled executable:
```bash
./real_time_eeg
```

The application provides:
- EEG signal visualization and processing
- TMS trigger control
- MRI/fMRI data visualization with ROI analysis
- Real-time phase estimation
- Data export capabilities

## Project Structure

- `UI/`: Qt-based user interface components
- `devices/`: Hardware interface implementations
  - `EEG/`: EEG device communication and processing
  - `TMS/`: TMS device control (MagPro)
- `dataHandler/`: Data processing and management
- `math/`: Mathematical operations and DSP functions
- `utils/`: Utility functions and helpers

## License

BSD 3-Clause License

Copyright (c) 2024 Veikka Leppänen, Copyright (c) 2024 Joonas Laurinoja, Copyright (c) 2024 Dogu Baran Aydogan.

See the [LICENSE](LICENSE.md) file for details.

## Performance Considerations

The application implements several optimizations:
- OpenMP parallel processing for computationally intensive operations
- Memory locking for real-time performance
- SIMD instructions (FMA) when available
- Efficient matrix operations using Eigen library

## Authors

- Veikka Leppänen
- Joonas Laurinoja
- Dogu Baran Aydogan