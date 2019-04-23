//#define DEBUGPRINT
/*#define MPI*/
/*#define OMP*/

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef MPI
#include <mpi.h>
#endif

#ifdef OMP
#include <omp.h>
#endif

//public functions
//return one pseudorandom number
double randomRSA();  

//return a vector of nrandomRSA double precision pseudorandom numbers 
int vrandomRSA(double * randomRSA, long nrandomRSA); 

//initialize state using time in microseconds since 1/1/1970, mpirank, and mpisize
int randomRSA_init();

//initialize state using seed, mpirank=0 and mpisize=1
int randomRSA_init_seed(uint64_t seed);

//initialize state using seed, mpirank, and mpisize
int randomRSA_init_seed_MPI(uint64_t seed, int mpirank, int mpisize);

//validate the operation of the vectorized pseudorandom number generator
//int randomRSA_validate(int numvectorcalls);

//get the current value of the exponent
int randomRSA_get_exponent(uint32_t *exponent);

//set exponent to odd value <=257. Recommend e<=17. If illegal, make no changes.
int randomRSA_set_exponent(uint32_t exponent);

//get the current values of the safe primes, as well as composite n=p1*p2
int randomRSA_get_primes(uint32_t *prime1, uint32_t *prime2, uint64_t *composite);

//set safe primes 2^32>prime1>prime2>2^31, and composite n=p1*p2. If illegal, make no changes.
int randomRSA_set_primes(uint32_t prime1, uint32_t prime2);

//get q=2^63-25, prime used in skip generator. q can not be changed)
int randomRSA_get_skipprime(uint64_t *skipprime);

//return the current value of the primitive root (mod q) (q=2^63-25 is prime used in skip generator and can not be changed)
int randomRSA_get_primitiveroot(uint32_t *primitiveroot);

//RSA vector size
int randomRSA_vectorsize();

//get total number of random numbers returned
long randomRSA_randsreturned();


//set primitive root (mod q), must be one of allowed values. If illegal, make no changes.
int randomRSA_set_primitiveroot(uint32_t primitiveroot);

int randomRSA_get_state(uint32_t *prime1, uint32_t *prime2, uint32_t *exponent, uint32_t *primitiveroot, uint32_t *index, uint32_t *hash, uint64_t *messages,uint64_t *skips);

int randomRSA_set_state(uint32_t prime1, uint32_t prime2, uint32_t exponent, uint32_t primitiveroot, uint32_t index, uint32_t hash, uint64_t *messages, uint64_t *skips);






//number of elements in vectors used in vectorized calculation 
//medium size value as tradeoff between memory use and speed
#define RSAVECTORSIZE 16384

//other possible values
//larger values give faster vector operations in exchange for more memory use 
//smaller values give slower vector operations in exchange for less memory use 
//#define RSAVECTORSIZE 1048576
//#define RSAVECTORSIZE 65536
//#define RSAVECTORSIZE 32768
//#define RSAVECTORSIZE 16384
//#define RSAVECTORSIZE 4096


//default exponent, larger than minimum value e=3, but small for speed
#define DEFAULTEXPONENT 5


