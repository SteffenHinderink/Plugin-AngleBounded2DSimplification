# Angle-Bounded 2D Mesh Simplification

This is the implementation for the paper "Angle-Bounded 2D Mesh Simplification" submitted to GMP 2022.
It is implemented as a plugin for the geometry processing software OpenFlipper.
It was tested on Mac and Linux.

## Prerequisites

- OpenFlipper: https://www.graphics.rwth-aachen.de/software/openflipper/
    - Tested with version 5.0

- Triangle: https://www.cs.cmu.edu/~quake/triangle.html
    - Tested with version 1.6

- Gurobi Optimizer: https://www.gurobi.com
    - Tested with version 9.5
    - Installation guide for Mac:
      https://www.gurobi.com/documentation/9.5/quickstart_mac/software_installation_guid.html#section:Installation
    - Installation guide for Linux:
      https://www.gurobi.com/documentation/9.5/quickstart_linux/software_installation_guid.html#section:Installation

## Build

- Move `triangle.h` and `triangle.c` from the Triangle library
  into this directory (``Plugin-AngleBounded2DSimplification``)

- Move this directory (``Plugin-AngleBounded2DSimplification``)
  into the OpenFlipper root directory (e.g. ``OpenFlipper-Free``)

- Build OpenFlipper in the OpenFlipper root directory (e.g. ``OpenFlipper-Free``):
    - ``mkdir build``
    - ``cd build``
    - ``cmake ..``
    - ``make``

- Run OpenFlipper:
    - Mac: ``./Build/OpenFlipper.app/Content/MacOS/OpenFlipper``
    - Linux: ``./Build/bin/OpenFlipper``

The plugin can be used through an interface in the toolbox on the right side of OpenFlipper.
Parameters for the triangulation and the decimation can be set.
