# pinverse

This repository contains a Fortran module for calculating the singular value decomposition (SVD) and pseudoinverse of a matrix.
============


## Table of Contents

- [pinverse](#pinverse)
- [This repository contains a Fortran module for calculating the singular value decomposition (SVD) and pseudoinverse of a matrix.](#this-repository-contains-a-fortran-module-for-calculating-the-singular-value-decomposition-svd-and-pseudoinverse-of-a-matrix)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [fpm](#fpm)
  - [Module Description](#module-description)
  - [Tests](#tests)
  - [Contributing](#contributing)
-----
## Requirements
To use the `pinverse` module, you need the following:

- Fortran compiler (Intel Compiler or gfortran)
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
- `pinverse`: Function to compute the pseudoinverse of a matrix using the SVD.
-----

## Tests

The tests directory contains test programs to verify the functionality of the pinverse module. To run the tests using fpm, you can use response files for specific compilers:

```bash
fpm @ifort
fpm @ifx
fpm @gfortran
```
-----

## Contributing

Contributions to pinverse are welcome! If you find any issues or would like to suggest improvements, please open an issue or submit a pull request.