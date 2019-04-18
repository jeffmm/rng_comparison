# A comparison of random number generators

This repo produces two binaries that produce random numbers given an initial seed. The purpose of the repo is to compare the header-only KISS RNG (based on UCL Professor David Jones' [JKISS implementation](http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf)) and the GSL RNG library. Comparisons were done in python and bash comparing random number autocorrelations, fitting to expected distributions, and speed.

Long story short, both generators properly returned pseudo-random numbers that passed my statistical tests, but KISS was better than two orders of magnitude faster than GSL, while having the advantage of being header-only.



