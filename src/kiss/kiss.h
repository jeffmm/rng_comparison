#ifndef _KISS_H_
#define _KISS_H_

#include <iostream>
#include <math.h>
#include <float.h>

/** KISS random number generator class. 
  * Based on UCL Professor David Jones' freely available JKISS RNG found at
  * http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf. Proper seeding
  * requires the initialization of four unsigned integer variables which are
  * used to generate subsequent random numbers.
  */
class KISSRNG {
  private:
    bool generate = false;
    double z1;
    unsigned int x = 123456789,
                 y = 987654321,
                 z = 43219876,
                 c = 6543217; 
  public:
    /** Initialize RNG using a single long seed. Two of the generator's four
     * variables are initialized by splitting the provided seed, and the
     * remaining two are provided for by the generator. This is an easier
     * method of initialization, but lacks the robustness of the main Init()
     * function. Since the variables initially have low entropy, this
     * initialization requires a "warm up", and immediately generates and
     * discards 1000 random numbers (see David Jones' discussion on warming up
     * RNGs).
     */
    void InitCold(long seed) {
      int seed1 = (int)(seed);
      int seed2 = (int)(seed >> 32);
      x = (seed1 < 0 ? -seed1 : seed1);
      z = (seed2 < 0 ? -seed2 : seed2);
      y = JKISS();
      c = JKISS();
      for (int i=0; i<1000; ++i) {
        JKISS();
      }
    }
    /** Initialize all four RNG variables with four unique unsigned integers.
     * There is a shorter "warm up" period of this generator since providing
     * four random unsigned integers is expected to have higher entropy. */
    void Init(unsigned int xx, unsigned int yy, 
              unsigned int zz, unsigned int cc) {
      x = xx;
      y = yy;
      z = zz;
      c = cc;
      for (int i=0; i<50; ++i) {
        JKISS();
      }
    }
    /** Generate a random unsigned integer. */
    unsigned int JKISS() {
      unsigned long long t;
      x = 314527869 * x + 1234567;
      y ^= y << 5; y ^= y >> 7; y ^= y << 22;
      t = 4294584393ULL * z + c; c = t >> 32; z = t;
      return x + y + z;
    }
    /** Returns random double in [0,1) with single precision. This method does
     * not have enough precision for sampling all possible doubles in the range
     * of [0,1), but more than good enough for single precision. If more
     * precision is needed, use the RandomUniformDbl method.
     */
    double RandomUniform() {
      double r = JKISS()/(UINT_MAX+1.0);
      return r;
    }
    /** Higher precision 53 bit random double in [0,1). Samples all possible
     * doubles in range, but requires a few more operations. */
    double RandomUniformDbl() {
      unsigned int a = JKISS() >> 6; /* Upper 26 bits */
      unsigned int b = JKISS() >> 5; /* Upper 27 bits */
      double r = (a * 134217728.0 + b) / 9007199254740992.0;
      return r;
    }
    /** Returns a random integer in the range [0,9]. Algorithm generates a
     * random 32 bit and returns the last integer in the sequence. Guaranteed
     * to be well-behaved for the JKISS generator.
    */
    int RandInt() {
      int ri = JKISS() % 10;
      return ri;
    }
    /** Generate a random variate from a normal distribution. Takes a mean mu
     * and standard deviation sigma as inputs. Algorithm uses the Box-Muller
     * transform. 
     */
    double RandomNormal(double mu, double sigma) {
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
    /** Generates a random unit vector in 2 or 3 dimensions. Takes the number
     * of dimensions and the vector as input.
     */
    void RandomUnitVector(int n_dim, double vect[]) {
        double x, y, z, w, t;

        if (n_dim < 2 || n_dim > 3) {
          std::cerr << "RandomUnitVector method requires n_dim = 2 or 3.\n";
          return;
        }
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
    /** Generates a vector of n_dim random numbers in the range [-0.5,0.5).
     * Takes the number of dimensions and vector as inputs. 
    */
    void RandomUniformVector(int n_dim, double vect[]) {
      for (int i=0; i<n_dim; ++i) {
        vect[i] = RandomUniform() - 0.5;
      }
    }
};

#endif // _KISS_H_

