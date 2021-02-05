// run in wsl, windows doesn't support gmp
#include <iostream>
#include <gmp.h>

// NPLEVELS, NPARAMS (from paper)
const double params[10][5] = {
    {80, 1600, 50, 1, 1.4}, {100, 5000, 70, 1, 1.45},
    {120, 8000, 90, 1, 1.5}, {140, 18000, 120, 1, 1.5},
    {160, 36000, 120, 1, 1.45}, {180, 66000, 120, 1, 1.45},
    {200, 120000, 150, 2 , 1.5}, {220, 200000, 180, 2 , 1.55},
    {230, 280000, 195, 2 , 1.55}, {240, 360000, 210, 2 , 1.57}
};

struct Sieve
{
    // number to factor
    mpz_t N;

    // upper bound for primes in factor base
    uint32_t fb_bound;

    // large prime bound (only relations whose large prime cofactors are smaller are added to hashtable)
    uint32_t lp_bound;

    // the sieve length
    uint32_t M;

    // the sieve threshold
    float T;

    Sieve(mpz_t n)
    {
        mpz_init_set (N, n);
        
        // select parameters
        int bits = mpz_sizeinbase (N, 2);

        std::cout << "choosing parameters for " << bits << " bit number... " << std::endl;
        
        int p1 = 0;
        int p2 = 0;
        double fac = 0;
        fb_bound = static_cast<uint32_t> (params[p1][1] * fac + params[p2][1] * (1 - fac));
        lp_bound = fb_bound * static_cast<uint32_t> (params[p1][2] * fac + params[p2][2] * (1 - fac));
        M = static_cast<uint32_t> (params[p1][3] * fac + params[p2][3] * (1 - fac));
        // todo add T
    }
};

int main(int argc, char *argv[])
{
    std::cout << "starting sieve..." << std::endl;
    
    mpz_t n;
	mpz_init (n);
    mpz_set_str (n, argv[1], 10);
    Sieve sieve(n);

    return 0;
}