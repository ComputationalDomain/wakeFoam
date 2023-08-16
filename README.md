# wakeFoam
An external library for OpenFOAM, containing a custom boundary condition for a velocity field behind cylinder.

<img src="https://img.shields.io/github/license/prabhuomkar/pytorch-cpp">

![OpenFOAM v2106](https://img.shields.io/badge/OpenFOAM-v2106-brightgreen.svg)
![OpenFOAM 8](https://img.shields.io/badge/OpenFOAM-6-brightgreen.svg)

## Table of Contents

1. [Introduction](#introduction)
2. [Implementaion](#implementation)
3. [License](#license)

## Introduction

WakeFoam library contains a custom inflow boundary condition, which uses a trained Neural Network to simulate the velocity field of a cylinder wake at the inflow boundary.

## Implementation

When using the wakeFoam library the ```net.pt``` and ```wake_BC``` files must be present in the simulations folder.

### Installation
```bash
cd $WM_PROJECT_USER_DIR
git clone https://github.com/ComputationalDomain/wakeFoam.git
cd wakeFoam
wmake
```

### Parameters
1. `flowSpeed` Defines the flowspeed outside of the wake (Note: the velocity should be adjusted according to parameter `D` such that the Reynolds number stays equal to 100).
2. `locationY` Defines the location of the center axis in spanwise direction.
3. `locationStreamwise` Defines the distance of inlet surface from the cylinder (e.g. 8 diameters behind cylinder).
4. `D` Diameter of the cylinder.
5. `streamwise` streamwise axis (e.g. x).
6. `spanwise` spanwise axis (e.g. y).

### Usage
There are tutorials located in `wakeFoam/tutorials`.

## License
This repository is licensed under MIT as given in [LICENSE](LICENSE).
