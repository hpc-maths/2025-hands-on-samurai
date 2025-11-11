:::{tip} How to start
If you have not already done so, activate the conda environment created in the [Setup Instructions](00-setup.md):

```bash
conda activate samurai-practical-session
```

For each step, we give you two main files, copy these two files into a new folder:
- `practical_session/material/step_name/CMakeLists.txt`
- `practical_session/material/step_name/main.cpp`

Then, create a build folder and compile the code as follows:

```bash
cmake -S . -B build
cmake --build build
cd build
./samurai_setup_test
```

We recommend changing the executable name in the `CMakeLists.txt` file to avoid confusion with the setup test and to keep track of your work at each step.
:::