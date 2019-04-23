#define N_RAND 1000000000
#include "randomRSA.h"

/** This program generates N_RAND random uniform numbers using the KISS RNG.
 * Expects a seed as an argument, otherwise runs with a default seed.
*/
int main(int argc, char * argv[]) {
  long seed;
  long i;
  double k;

  seed = 77777777777777;
  /*if (argc == 1) {*/
    /*printf("Running with default seed: %ld\n", seed);*/
  /*}*/
  /*else {*/
    /*seed = (long) argv[1];*/
  /*}*/
   /*Initialize the RNG using the seed*/
  k = 0;
  randomRSA_init();
  
  for (i=0; i<N_RAND; ++i) {
    k = randomRSA();
     /*For generating the numbers to file. Ruins time benchmarks.*/
  }
  printf("%2.8f\n", k);

  return 0;
}
