# nBody and Orbital Simulation

This repository features two main files:

- **nbody.cpp**: A C++ implementation of an N-body simulation with both O(N²) direct summation and faster O(N log N) methods (e.g., Barnes-Hut tree code).
- **ca_ode_orbits.py**: A Python script using various integrators (Explicit Euler, Runge-Kutta, Leap Frog, etc.) to simulate orbital paths around a central mass.

## Requirements

- C++17 compiler
- Python 3 (for the orbital simulation)
- (Optional) CMake or Make for building the C++ project
- (Optional) Additional plotting libraries such as matplotlib (Python)

## Building nBody

```bash
cmake -B build
cmake --build build
```

Or simply:

```bash
make
```

This creates the `nbody.exe` executable on Windows.

## Running nBody

Run through the terminal:

```bash
./nbody.exe
```

The advanced force calculations (Barnes-Hut or alternative O(N log N) methods) can be toggled via code-level flags/functions.

## Using ca_ode_orbits.py

Run the Python script for orbital integrator comparisons:

```bash
python ca_ode_orbits.py
```

It demonstrates multiple integration schemes and plots the results with matplotlib.

## Further Insights

- **nbody.cpp**:
  - Direct summation O(N²) and tree-based O(N log N) force approximations.
  - Additional routines for energy calculation, half-mass radius, and shell-based analysis.
- **ca_ode_orbits.py**:
  - Integrators (Explicit Euler, Runge-Kutta 2/4, Leap Frog, Semi-Implicit Euler).
  - Compares orbital accuracy over several revolutions.
