#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "../src/rng/rng.h"


TEST_CASE("GSL RNG functions") {
  RNG rng;
  rng.Init(rand());
  SECTION("Random number ranges are defined") {
    REQUIRE(rng.RandomUniform() < 1);
    REQUIRE(rng.RandomUniform() >= 0);
  }
}

