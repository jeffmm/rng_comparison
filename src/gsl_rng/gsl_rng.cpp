#define N_RAND 1000000000
#include "gsl_rng.h"

/** This program generates N_RAND random uniform numbers using the GSL RNG.
 * Expects a seed as an argument, otherwise runs with a default seed.
*/
int main(int argc, char * argv[]) {
  long seed = 77777777777777;
  if (argc == 1) {
    printf("Running with default seed: %ld\n", seed);
  }
  else {
    try {
      seed = (long) argv[1];
    }
    catch (...) {
      // In case we can't format the argument to long
      printf("ERROR: Could not parse seed argument\n");
      exit(1);
    }
  }
  // Initialize the RNG using the seed
  GSLRNG rng(seed);
  double k = 0;

  for (int i=0; i<N_RAND; ++i) {
    k = rng.RandomUniform();
    // For generating the numbers to file. Ruins time benchmarks.
    //std::cout << k <<"\n";
  }
  // Keep the compiler from skipping the loop
  std::cout << k << "\n";

  return 0;
}
