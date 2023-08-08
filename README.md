[![GitHub](https://img.shields.io/badge/GitHub-ForInv-blue.svg?style=social&logo=github)](https://github.com/gha3mi/forinv)
[![Version](https://img.shields.io/github/v/tag/gha3mi/forinv?color=blue&logo=github&style=flat)](https://github.com/gha3mi/forinv/releases)
[![Documentation](https://img.shields.io/badge/ford-Documentation%20-blueviolet.svg)](https://gha3mi.github.io/forinv/)
[![License](https://img.shields.io/github/license/gha3mi/forinv?color=green)](https://github.com/gha3mi/forinv/blob/main/LICENSE)
[![Build](https://github.com/gha3mi/forinv/actions/workflows/ci.yml/badge.svg)](https://github.com/gha3mi/forinv/actions/workflows/ci.yml)

<img alt="ForInv" src="https://github.com/gha3mi/forinv/raw/main/media/logo.png" width="750">

**ForInv**: A Fortran library for inverse and pseudo-inverse calculations.


## Requirements
To use the `forinv` module, you need the following:

- Fortran compiler
- LAPACK or MKL

## fpm dependency

If you want to use `ForInv` as a dependency in your own fpm project,
you can easily include it by adding the following line to your `fpm.toml` file:

```toml
[dependencies]
forinv = {git="https://github.com/gha3mi/forinv.git"}
```

## Usage

Here is an example of how to use the `forinv` module in your Fortran code:
```fortran
program main

   use kinds
   use forinv, only: inv

   implicit none

   ! Declare variables
   real(rk), dimension(:, :), allocatable :: A, A_inv

   ! Initialize matrix A
   ...

   ! Call pseudoinverse function
   A_inv = inv(A)

end program main
```

## How to run tests and examples

**Clone the repository:**

You can clone the `ForInv` repository from GitHub using the following command:

```shell
git clone https://github.com/gha3mi/forinv.git
```

```shell
cd forinv
```

**Run tests:**

To set the stack size to unlimited, use the following command: `ulimit -s unlimited`.

```shell
fpm @gfortran-test
```

```shell
fpm @ifort-test
```

```shell
fpm @ifx-test
```

```shell
fpm @nvfortran-test
```

## API documentation

The most up-to-date API documentation for the master branch is available
[here](https://gha3mi.github.io/forinv/).
To generate the API documentation for `ForInv` using
[ford](https://github.com/Fortran-FOSS-Programmers/ford) run the following
command:

```shell
ford ford.yml
```

## Contributing
Contributions to `ForInv` are welcome! If you find any issues or would like to suggest improvements, please open an issue.
