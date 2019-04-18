#ifndef _RNG_GSL_H_
#define _RNG_GSL_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <iostream>

/** Random number generator class using GSL libraries. Memory allocation and
 * cleanup functions are handled by the class.
  */
class GSLRNG { 
  private:
    const gsl_rng_type *T;
    void Clear() {
      gsl_rng_free(r);
    }
    gsl_rng *r;
  public:
    GSLRNG() {}
    /** Instantiation of RNG with a default seed value. Calls Init(seed).
    */
    GSLRNG(long seed) {
      Init(seed);
    }
    ~GSLRNG() {
      Clear();
    }
    /** Initialize the RNG with a default seed. Accepts a long for the seed.
      */
    void Init(long seed) {
      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc(T);
      gsl_rng_set(r, seed);
    }
    GSLRNG(const GSLRNG& that) {
      this->Init(gsl_rng_get(that.r));
    }
    GSLRNG& operator=(GSLRNG const&that) {
      this->Init(gsl_rng_get(that.r));
      return *this;
    }
    /** Returns a random uniform in the range [0,1). Uses the GSL function
     * gsl_rng_uniform on the generator.
      */
    double RandomUniform() {
      return gsl_rng_uniform(r);
    }
};

#endif // _RNG_GSL_H_
