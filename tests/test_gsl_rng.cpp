#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "../src/gsl_rng/gsl_rng.h"


TEST_CASE("GSL RNG functions") {
  GSLRNG rng;
  rng.Init(rand());
  SECTION("Random number ranges are defined") {
    REQUIRE(rng.RandomUniform() < 1);
    REQUIRE(rng.RandomUniform() >= 0);
  }
}

