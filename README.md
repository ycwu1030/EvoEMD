# EvoEMD : cosmic Evolution with an Early Matter-Dominated era

A framework to evaluate the evolution w/ or w/o an early matter dominated (EMD) era.

## Features

- Global parameter system
- Global particle system
- Global process system
- Automatically build up the Boltzmann equation according to the user's definition of particle and process
- Solving the Boltzmann equation using 4th order Runge-Kutta method with adaptive steps tailored to Cosmology application
- Cache the collision rate calculation results for fast evaluation

## Document

In preparation...

## Dependencies

- [GSL](https://www.gnu.org/software/gsl/)
- [Cuba](http://www.feynarts.de/cuba/)

## Installation

Make sure that `GSL` is there which can be installed through many package management systems (`apt-get`, `brew`, etc.). Downloading `Cuba` from the website and compiling it such that one can find `cuba.h` and `libcuba.a` in its directory.

```bash
git clone https://github.com/ycwu1030/EvoEMD.git
cd EvoEMD
mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=INSTALL_PATH -DCUBA_ROOT=CUBA_PATH ../
make
make install # Optional
```
where `INSTALL_PATH` should be replaced by the directory where you want to install the `EVOEMD` library, `CUBA_PATH` should be replaced by the `CUBA` directory. Then one can find the executable for the examples in `EvoEMD/build/bin`.

The last command is optional, which install the whole library is into `INSTALL_PATH` (If it is a system folder, one should also use `sudo`). The headers are in `INSTALL_PATH/include`, while the shared library is in `INSTALL_PATH/lib`. If users want to integrate `EvoEMD` into their own CMake project, we also provide `FindEVOEMD.cmake` in `INSTALL_PATH/share/cmake/Modules` as well as `SOURCE_DIR/cmake/Modules`.

### Examples

Simple examples can be found in corresponding folders under `Models`.
