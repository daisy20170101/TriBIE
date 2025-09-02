#!/bin/bash

echo "=========================================="
echo "Compiling TriGreen (local DP version)..."
echo "=========================================="

# Compile the main program directly (no separate modules needed)
echo "Compiling calc_trigreen.f90..."
mpif90 -o calc_trigreen calc_trigreen.f90

if [ $? -eq 0 ]; then
    echo "✓ Compilation successful!"
    echo "Executable: calc_trigreen"
else
    echo "✗ Compilation failed!"
    exit 1
fi

echo "=========================================="
echo "Compilation completed!"
echo "=========================================="
