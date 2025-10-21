# Setup Instructions

During this practical session, we will be using the samurai library and demonstrating its capabilities through a series of coding exercises. Please follow the instructions below to set up your development environment.

## Prerequisites

We provide a conda environment file `environment.yml` to facilitate the setup process. You can create and activate the conda environment by running the following commands in your terminal:

```bash
conda env create -f conda/environment.yml
conda activate samurai-practical-session
```

If you want to set up the environment manually, please ensure you have the following dependencies installed:

- A C++ compiler that supports C++20 (e.g., GCC, Clang)
- CMake (version 3.16 or later)
- xtensor library (version 0.26 or later)
- HighFive library (version 3.0 or later)
- MPI library (e.g., OpenMPI or MPICH)
- fmt library (version 11.0 or later)
- pugixml library (version 1.15 or later)
- cli11 library (version less than 2.5)

## Verifying the Setup

To verify that your environment is set up correctly, you can compile and run the provided example code. Navigate to the `practical_session/material/00-setup` directory and execute the following commands:

```bash
cmake -S . -B build
cmake --build build
./build/samurai_setup_test -h
```

You should see the help message for the samurai setup test program, indicating that everything is working correctly.

If you encounter any issues during the setup process, please refer to the [samurai documentation](https://hpc-math-samurai.readthedocs.io/) or reach out to the course instructors for assistance.

You are now ready to begin the practical session! Enjoy coding with samurai!