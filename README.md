# A comparison of random number generators

This repo produces two binaries that produce random numbers given an initial seed. The purpose of the repo is to compare the header-only KISS RNG (based on UCL Professor David Jones' [JKISS implementation](http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf)) and the GSL RNG library. Comparisons were done in python and bash comparing random number autocorrelations, fitting to expected distributions, and speed.

Long story short, both generators properly returned pseudo-random numbers that passed my statistical tests, but KISS was over 4 times faster than GSL, while having the advantage of being header-only.

## Compiling

Requirements are cmake, a C++11 compiler, and the GSL library installed somewhere where cmake can find it (like */usr/local/lib*, etc).

To compile using cmake, do the usual 

```
mkdir build
cd build
cmake ..
make
```

## Documentation

By default, documentation will be built into the *build/documentation* folder in both html and latex using Doxygen (if it is available).

## Testing

To run unit tests, navigate to a fresh *build* folder and build the binaries in debug mode with `cmake -DCMAKE_BUILD_TYPE=Debug ..`, and the unit test binaries will build as well. Running `ctest` or `make test` will run all unit test binaries and report any failures. If any of the tests fail, the offending unit test binary should be run with `./tests/test_<rng_name>` to see which exact unit test failed, using the header-only Catch2 library.

