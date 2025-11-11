# Environment Setup

:::{note} Main Objectives
- Set up the development environment for the practical session
- Verify that samurai and visualization tools are working correctly
:::

During this practical session, we will use the samurai library and demonstrate its capabilities through a series of coding exercises. Please follow the instructions below to set up your development environment.

## Prerequisites

Before starting, ensure you have the following installed on your system:
- **Git** (to clone the repository)
- **Conda**, **Miniconda** or **Micromamba** (recommended for easy setup)
- **ParaView** (for visualization of 2D/3D results)

If you prefer a manual setup without conda, you will also need:
- A C++ compiler supporting C++20 (e.g., GCC 10+, Clang 12+)
- CMake (version 3.16 or later)
- The following libraries: xtensor, HighFive, MPI, fmt, pugixml, cli11 (see details below)

The first step is to clone the samurai hands-on repository from GitHub. You can do this by running the following command in your terminal:

```bash
git clone https://github.com/hpc-maths/2025-hands-on-samurai.git
```

Go to the cloned directory:

```bash
cd 2025-hands-on-samurai
```

We provide a conda environment file (`environment.yml`) to facilitate the setup process. You can create and activate the conda environment by running the following commands in your terminal:

```bash
conda env create -f conda/environment.yml
conda activate samurai-practical-session
```

If you want to set up the environment manually, please ensure you have the following dependencies installed:

| Dependency   | Minimum Version | Notes                             |
| ------------ | --------------- | --------------------------------- |
| C++ compiler | C++20 support   | GCC 10+, Clang 12+, or MSVC 2019+ |
| CMake        | 3.16            |                                   |
| xtensor      | 0.26            |                                   |
| HighFive     | 3.0             | HDF5 C++ interface                |
| MPI          | -               | OpenMPI or MPICH                  |
| fmt          | 11.0            | String formatting library         |
| pugixml      | 1.15            | XML parsing library               |
| cli11        | 2.4.x           | Version 2.5+ not compatible       |

## Verifying the setup

To verify that your environment is set up correctly, you can compile and run the provided example code. Navigate to the `practical_session/material/00-setup` directory and execute the following commands:

```bash
cd practical_session/material/00-setup
cmake -S . -B build
cmake --build build
./build/samurai_setup_test -h
```

**Expected output:** You should see a help message listing the available command-line options for the test program. If you see this, your samurai installation is working correctly!

:::{caution}Common issues
- cmake not found: Make sure cmake is in your PATH or activate the conda environment
- Compiler errors: Ensure your compiler supports C++20
- Missing libraries: Verify all dependencies are installed (run conda list if using conda)
:::

If you encounter any issues during the setup process, please refer to the [samurai documentation](https://hpc-math-samurai.readthedocs.io/) or reach out to the course instructors for assistance.

:::{important}
After opening a new terminal later, remember to run:

```bash
conda activate samurai-practical-session
```
:::

## Verify the visualization tools

### For Multi-Dimensional Simulations (2D/3D)

To visualize results from 2D or 3D simulations, you will need ParaView. Download it from the [official ParaView website](https://www.paraview.org/download/) if not already installed.

**Testing ParaView visualization:**

1. Open ParaView
2. File → Open → Navigate to `practical_session/material/00-setup/outputs/2d_example.xdmf`
3. Select the `XDMF Reader` in the list
4. Click **Apply** in the Properties panel
5. You should see a 2D field visualization

:::{tip}
Always open the `.xdmf` file (not the `.h5` file) to load both mesh and field data.
:::

### For 1D simulations

ParaView does not support 1D visualizations directly. Instead, use the Python script provided with samurai.

**Testing 1D visualization:**

Navigate to the outputs directory and run the script:

```bash
cd practical_session/material/00-setup/outputs
python read_mesh.py 1d_example --field u
```

**Expected output:** A matplotlib window should appear showing a 1D plot of the field `u`.

:::{caution}Common issues
- matplotlib not installed: Run conda install matplotlib or pip install matplotlib
- No display: If running on a remote server, use `--save filename` option to save the plot instead (the png extension is automatically added to the filename)
:::

## Next Steps

If you successfully:
- ✅ Built and ran the test program
- ✅ Visualized the 2D example in ParaView
- ✅ Visualized the 1D example with the Python script

**Congratulations!** 🎉 Your environment is fully set up and you're ready to proceed with the practical sessions.

**Troubleshooting:** If you encountered any issues, please:
- Check the [samurai documentation](https://hpc-math-samurai.readthedocs.io/)
- Review the error messages carefully
- Ask the course instructors for assistance

Now, proceed to [First Steps](01-first-step.md) to begin the practical session. Enjoy coding with samurai!
