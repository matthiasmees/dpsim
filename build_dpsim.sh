#!/bin/bash

# Step 1: Build file
cd /dpsim/build
# cmake --build . --target EMT_cosim_9bus_4order -j4
# cmake --build . --target DP_cosim_9bus -j4
# cmake --build . --target cosim-9bus -j4
cmake --build . --target DP_RXLoad -j4


