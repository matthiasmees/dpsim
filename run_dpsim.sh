#!/bin/bash

# Step 1: Build file
cd /dpsim/build
cmake --build . --target EMT_cosim_9bus_4order -j4
cmake --build . --target DP_cosim_9bus_4order -j4

# Step 2: Export Python path
cd /dpsim/build
export PYTHONPATH=$(pwd):$(pwd)/../python/src/

# Step 3: Open Python Notebooks
cd /dpsim
jupyter lab --ip="0.0.0.0" --allow-root --no-browser
