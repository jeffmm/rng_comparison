#define N_RAND 100000
#include "kiss.h"

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
      printf("ERROR: Could not parse seed argument\n");
      exit(1);
    }
  }
  RNG rng;
  rng.InitCold(seed);
  //double rand[N_RAND];
  //double k = 0;
  
  for (int i=0; i<N_RAND; ++i) {
    //k = rng.RandomUniform();
    std::cout << rng.RandomUniform() <<"\n";
  }

  return 0;

}
