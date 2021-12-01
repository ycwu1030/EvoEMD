# EvoEMD : cosmic Evolution with an Early Matter-Dominated era

A framework to evaluate the evolution w/ or w/o an early matter dominated (EMD) era.

## Features

- Global parameter system
- Global particle system
- Global process system
- Different methods for Hubble parameter calculation
- Automatically build up the Boltzmann equation according to the user's definition of particle and process
- Solving the Boltzmann equation using 4th order Runge-Kutta method with adaptive steps tailored to Cosmology application
- Cache the collision rate calculation results for fast evaluation

## Document

The physics behind the `EvoEMD` and technical details can be found in:

- M. Dutra, Y. Wu; EvoEMD: cosmic evolution with an early matter-dominated era; [arXiv: 2111.15665](https://arxiv.org/abs/2111.15665)

If you are using `EvoEMD`, please cite above paper.

## Dependencies

- [GSL](https://www.gnu.org/software/gsl/) (At least v2.0 is required)
- [Cuba](http://www.feynarts.de/cuba/)
- [spdlog](https://github.com/gabime/spdlog) (Included in `EvoEMD`)

## Installation

Make sure that `GSL` (later than v2.0) is there which can be installed through many package management systems (`apt-get`, `brew`, etc.). Downloading `Cuba` from the website and compiling it such that one can find `cuba.h` and `libcuba.a` in its directory.

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

Simple examples can be found in corresponding folders under `Models`:
- ToyDM
    - A toy DM model, including freeze-in and freeze-out case.
- ToyLeptogenesis
    - A toy Leptogenesis model, the heavy neutrino can be initially thermalized or un-thermalized.
