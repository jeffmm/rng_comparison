#ifndef _RNG_GSL_H_
#define _RNG_GSL_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <iostream>

class RNG { 
  private:
    const gsl_rng_type *T;
    void Clear() {
      gsl_rng_free(r);
    }
    gsl_rng *r;
  public:
    RNG() {}
    RNG(long seed) {
      Init(seed);
    }
    ~RNG() {
      Clear();
    }
    void Init(long seed) {
      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc(T);
      gsl_rng_set(r, seed);
    }
    RNG(const RNG& that) {
      this->Init(gsl_rng_get(that.r));
    }
    RNG& operator=(RNG const&that) {
      this->Init(gsl_rng_get(that.r));
      return *this;
    }
    double RandomUniform() {
      return gsl_rng_uniform_pos(r);
    }
};

#endif // _RNG_GSL_H_
