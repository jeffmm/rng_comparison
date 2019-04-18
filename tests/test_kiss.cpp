#define CATCH_CONFIG_MAIN
#include <catch.hpp>

/* Test RNG utilities to ensure that the unit vector is indeed a unit vector,
 * and the uniform vector consists of random numbers r such that -0.5 < r < 0.5
 * */
#include "../src/kiss/kiss.h"


TEST_CASE("KISS RNG functions") {
  KISSRNG rng;
  rng.Init(rand(), rand(), rand(), rand());
  double vec[3] = {0,0,0};
  SECTION("Random number ranges are defined") {
    REQUIRE(rng.JKISS() <= UINT_MAX);
    REQUIRE(rng.JKISS() > 0);
    REQUIRE(rng.RandomUniform() < 1);
    REQUIRE(rng.RandomUniform() >= 0);
    REQUIRE(rng.RandInt() < 10);
    REQUIRE(rng.RandInt() >= 0);
  }
  SECTION("RandomUnitVector returns unit vectors") {
    rng.RandomUnitVector(3, vec);
    REQUIRE(abs(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] - 1.0) < 1e-12);
    rng.RandomUnitVector(2, vec);
    REQUIRE(abs(vec[0]*vec[0] + vec[1]*vec[1] - 1.0) < 1e-12);
  }
  SECTION("RandomUniformVector returns vectors with entries r_i, such that -0.5 <= r_i < 0.5") {
    rng.RandomUniformVector(3, vec);
    REQUIRE(abs(vec[0]) <= 0.5);
    REQUIRE(abs(vec[1]) <= 0.5);
    REQUIRE(abs(vec[2]) <= 0.5);
  }
}

