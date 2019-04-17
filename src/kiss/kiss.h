#ifndef _KISS_H_
#define _KISS_H_

#include <iostream>
#include <math.h>
#include <float.h>

/* Public domain code for JKISS RNG found here:
 * http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf */

class RNG {
  private:
    bool generate = false;
    double z1;
    unsigned int x = 123456789,
                 y = 987654321,
                 z = 43219876,
                 c = 6543217; 
  public:
    /* Initialize from single uint and warm up the generator */
    void InitCold(long seed) {
      int seed1 = (int)(seed);
      int seed2 = (int)(seed >> 32);
      x = (seed1 < 0 ? -seed1 : seed1);
      z = (seed2 < 0 ? -seed2 : seed2);
      y = JKISS();
      c = JKISS();
      for (int i=0; i<100; ++i) {
        JKISS();
      }
    }
    /* Initialize with four uints */
    void Init(unsigned int xx, unsigned int yy, 
              unsigned int zz, unsigned int cc) {
      x = xx;
      y = yy;
      z = zz;
      c = cc;
    }
    unsigned int JKISS() {
      unsigned long long t;
      x = 314527869 * x + 1234567;
      y ^= y << 5; y ^= y >> 7; y ^= y << 22;
      t = 4294584393ULL * z + c; c = t >> 32; z = t;
      return x + y + z;
    }
    /* Returns random double in [0,1). Not enough precision for
     * sampling all possible doubles in that range, but more than
     * good enough for single precision.*/
    double RandomUniform() {
      double r = JKISS()/(UINT_MAX+1.0);
      return r;
    }
    /* Higher precision 53 bit random double in [0,1). Samples all possible
     * doubles in range.*/
    double RandomUniformDbl() {
      unsigned int a = JKISS() >> 6; /* Upper 26 bits */
      unsigned int b = JKISS() >> 5; /* Upper 27 bits */
      double r = (a * 134217728.0 + b) / 9007199254740992.0;
      return r;
    }
    /* Returns random integer in [0,9) */
    int RandInt() {
      int ri = JKISS() % 10;
      return ri;
    }
    double RandomGaussian(double mu, double sigma) {
      generate = !generate;
      if (!generate) {
         return z1 * sigma + mu;
      }
      double u1, u2;
      do {
         u1 = RandomUniform();
         u2 = RandomUniform();
       } while ( u1 <= DBL_MIN );

      double z0;
      z0 = sqrt(-2.0*log(u1)) * cos(2.0*M_PI * u2);
      z1 = sqrt(-2.0*log(u1)) * sin(2.0*M_PI * u2);
      return z0 * sigma + mu;
    }
    void RandomUnitVector(int n_dim, double vect[]) {
        double x, y, z, w, t;

        w = 1.0;
        if (n_dim == 3) {
            z = 2.0 * RandomUniform() - 1.0;
            w = sqrt(1 - z * z);
            vect[2] = z;
        }

        t = 2.0 * M_PI * RandomUniform();
        x = w * cos(t);
        y = w * sin(t);
        vect[0] = x;
        vect[1] = y;
    }
    void RandomUniformVector(int n_dim, double vect[]) {
      for (int i=0; i<n_dim; ++i) {
        vect[i] = RandomUniform() - 0.5;
      }
    }
};

#endif // _KISS_H_

