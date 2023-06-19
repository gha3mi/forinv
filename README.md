# pinverse
This repository contains a Fortran module for calculating the pseudoinverse of a matrix using singular value decomposition (SVD).

-----


## Table of Contents

- [pinverse](#pinverse)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [fpm](#fpm)
  - [Module Description](#module-description)
  - [Usage](#usage)
  - [Tests](#tests)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
-----
## Requirements
To use the `pinverse` module, you need the following:

- Fortran compiler (tested with Intel and NVIDIA)
- LAPACK or MKL

## Installation

### fpm
pinverse can be cloned and then built using [fpm](https://github.com/fortran-lang/fpm), following the instructions provided in the documentation available on Fortran Package Manager.

```bash
git clone https://github.com/gha3mi/pinverse.git
cd pinvers
fpm install --perfix .
```

Or you can easily include this package as a dependency in your `fpm.toml` file.

```toml
[dependencies]
[dependencies.pinverse]
git = "https://github.com/gha3mi/pinverse.git"
```

-----
## Module Description

The `pinverse` module provides functions and subroutines for calculating the SVD and pseudoinverse of a matrix. It includes the following functionalities:

- `svd`: Subroutine to compute the SVD of a matrix.
- `pinv`: Function to compute the pseudoinverse of a matrix using the SVD.
-----

## Usage
Here is an example of how to use the `pinverse` module in your Fortran code:
```fortran
program main

   use :: kinds
   use :: pinverse, only: pinv

   implicit none

   ! Declare variables
   real(rk), dimension(:, :), allocatable :: A, A_pinv

   ! Initialize matrix A
   ...

   ! Call pseudoinverse function
   A_pinv = pinv(A)

end program main
```
-----

## Tests

The `tests` directory contains test programs to verify the functionality of the `pinverse` module. To run the tests using `fpm`, you can use response files for specific compilers:

- For Intel Fortran Compiler (ifort):
```bash
fpm @ifort
```
Compiler flags: ```-Ofast -xHost -mtune=native -parallel -qmkl=parallel```

- For Intel Fortran Compiler (ifx):
```bash
fpm @ifx
```
Compiler flags: ```-Ofast -xHost -mtune=native -parallel -qmkl=parallel```

- For NVIDIA Compiler (nvfortran):
```bash
fpm @nvidia
```
Compiler flags: ```-O4 -mtune=native -llapack```

- For GNU Fortran Compiler (gfortran):
```bash
fpm @gfortran
```
Compiler flags ```-Wno-line-truncation -Ofast -march=native -llapack -lblas```

-----

## Documentation
To generate the documentation for the `pinverse` module using [ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following command:
```bash
ford project.yml
```

-----

## Contributing

Contributions to pinverse are welcome! If you find any issues or would like to suggest improvements, please open an issue or submit a pull request.