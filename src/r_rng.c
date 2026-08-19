/*
 * r_rng.c — R-compatible uniform random number stream.
 *
 * Reproduces R's default RNG exactly: set.seed(seed) followed by
 * unif_rand() draws from the Mersenne-Twister, including R's LCG seed
 * scrambling and its fixup of boundary values.  Needed so that seeded
 * papa2 runs emit the same bootstrap samples as R sessions that call
 * set.seed() before assignTaxonomy().
 */

#include <stdint.h>

#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL

typedef struct {
    uint32_t mt[MT_N];
    int mti;
} RRng;

/* R's set.seed: initial scrambling then per-word LCG fill (RNG.c: RNG_Init). */
void r_rng_seed(RRng *rng, uint32_t seed) {
    int j;
    for (j = 0; j < 50; j++) seed = 69069 * seed + 1;
    for (j = 0; j < MT_N + 1; j++) {
        seed = 69069 * seed + 1;
        if (j == 0) {
            /* i_seed[0] is the MT position word; FixupSeeds sets it to 624 */
            continue;
        }
        rng->mt[j - 1] = seed;
    }
    rng->mti = MT_N;
}

static uint32_t mt_genrand(RRng *rng) {
    uint32_t y;
    static const uint32_t mag01[2] = {0x0UL, MATRIX_A};
    uint32_t *mt = rng->mt;

    if (rng->mti >= MT_N) {
        int kk;
        for (kk = 0; kk < MT_N - MT_M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < MT_N - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        rng->mti = 0;
    }

    y = mt[rng->mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}

/* R's unif_rand for the Mersenne-Twister kind, including fixup(). */
double r_rng_unif(RRng *rng) {
    /* MT_genrand returns reals in [0,1) via the 2.3283064365386963e-10
     * multiplier; fixup() then keeps the value strictly inside (0,1). */
    double value = mt_genrand(rng) * 2.3283064365386963e-10;
    if (value <= 0.0) return 0.5 * 2.328306437080797e-10;
    if (1.0 - value <= 0.0) return 1.0 - 0.5 * 2.328306437080797e-10;
    return value;
}

/* Fill out[0..n) with a seeded R runif stream (helper for one-shot use). */
void r_rng_runif_fill(uint32_t seed, double *out, long long n) {
    RRng rng;
    long long i;
    r_rng_seed(&rng, seed);
    for (i = 0; i < n; i++) out[i] = r_rng_unif(&rng);
}
