# A comparison of random number generators

This repo produces two binaries that produce random numbers given an initial seed. The purpose of the repo is to compare the header-only KISS RNG (based on UCL Professor David Jones' [JKISS implementation](http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf)) and the GSL RNG library. Comparisons were done in python and bash comparing random number autocorrelations, fitting to expected distributions, and speed.

Long story short, both generators properly returned pseudo-random numbers that passed my statistical tests, but KISS was better than two orders of magnitude faster than GSL, while having the advantage of being header-only.

## Compiling

To compile using cmake, do the usual 

```
mkdir build
cd build
cmake ..
make
```

This will require a C++11 compiler and the GSL library installed somewhere where cmake can find it (like /usr/local/lib, etc).

## Testing

Running `ctest` or `make test` will run all unit tests (written using the header-only catch2 library) for each binary. If either of the tests fail, the unit test binaries themselves should be run with `./tests/test_<rng_name>` to see which unit test failed.

