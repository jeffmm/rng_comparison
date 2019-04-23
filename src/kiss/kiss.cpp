#define N_RAND 1000000000
#include "kiss.h"

/** This program generates N_RAND random uniform numbers using the KISS RNG.
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
  KISSRNG rng;
  rng.InitCold(seed);

  double k = 0;
  
  for (long i=0; i<N_RAND; ++i) {
    k = rng.RandomUniform();
    k+=1;
    // For generating the numbers to file. Ruins time benchmarks.
    //std::cout << k <<"\n"; 
  }
  printf("%2.2f\n",k);
  return 0;
}
