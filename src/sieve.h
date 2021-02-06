#pragma once

#include <algorithm>
#include <unordered_map>
#include <mutex>
#include "params.h"
#include <array>
#include <gmpxx.h>
#include <gmp.h>
#include <math.h>

const std::array<Params, 10> PARAMS_ARR = {
    Params(80,  1600, 50,  1, 1.4),
    Params(100,  5000,  70, 1, 1.45),
    Params(120,  8000,  90,  1 , 1.5),
    Params(140, 18000,  120, 1 , 1.5),
    Params(160, 36000,  120, 1 , 1.45),
    Params(180, 66000,  120, 1 , 1.45),
    Params(200, 120000, 150, 2 , 1.5),
    Params(220, 200000, 180, 2 , 1.55),
    Params(230, 280000, 195, 2 , 1.55),
    Params(240, 360000, 210, 2 , 1.57)
};

const std::array<uint32_t, 18> SMALL_PRIMES = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61
};

struct Sieve
{
    std::mutex lock;

    mpz_class &nc;
    uint32_t fb_bound;
    uint32_t lp_bound;
    uint32_t M;
    float T;
    uint32_t multiplier;
    std::vector<uint32_t> fb;
    std::vector<uint32_t> roots;
    std::vector<uint8_t> fb_logs;
    
    uint32_t ntheads = 1;
    uint32_t info_npoly = 0;
    uint32_t info_npg = 0;
    uint32_t nfull = 0;
    uint32_t npartial = 0;
    uint32_t tdiv_ct = 0;
    uint64_t sieve_locs = 0;
    uint32_t extra_rels = 120;

    double get_multiplier_score(uint32_t mult) const {
        mpz_class tempc;
        mpz_mul_ui(nc.get_mpz_t(), nc.get_mpz_t(), mult);
        double res = -0.5 * log(mult);
        if (mpz_mod_ui(tempc.get_mpz_t(), nc.get_mpz_t(), 8) == 1) {
            res += 2 * log(2.0);
        }

        for (size_t i = 1; i < 18; i++) {
            if (SMALL_PRIMES[i] == mult) {
                res += (1.0/SMALL_PRIMES[i]) * log(SMALL_PRIMES[i]);
            }
            else if (mpz_kronecker_ui(nc.get_mpz_t(), SMALL_PRIMES[i]) == 1) {
                res += (2.0/SMALL_PRIMES[i]) * log(SMALL_PRIMES[i]);
            }
        }

        mpz_divexact_ui(nc.get_mpz_t(), nc.get_mpz_t(), mult);
        return res;
    }

    std::vector<char> era_sieve(uint32_t fbb) const {
        std::vector<char> vals(fbb);
        for (int skip = 2; skip < static_cast<int>(sqrt(fbb) + 1); skip++) {
            if (vals[skip-2] == 1) {
                continue;
            }
            for (int pos = 2*skip; pos < fb_bound; pos += skip) {
                vals[pos-2] = 1;
            }
        }
        return std::move(vals);
    }

    std::vector<uint32_t> extract_primes(uint32_t fbb, uint32_t mult, const std::vector<char> &vals) const {
        std::vector<uint32_t> primes;
        for (int i = 0; i < fbb-2; i++) {
            if (vals[i] == 0 && (i == 0 || i+2 == mult || mpz_kronecker_ui(nc.get_mpz_t(), i+2) == 1)) {
                primes.push_back(i+2);
            }
        }
        return std::move(primes);
    }

    uint8_t fast_log (uint32_t x) const {
        x = (x*3)/2;
        uint8_t res = 0;
        while (x > 0){
            x = x >> 1;
            res++;
        }
        return res - 1;
    }

    uint32_t find_root_pocklington(uint32_t p, mpz_class &temp) {
        mpz_class pol;
        mpz_class temp2;
	    if (p % 4 == 3) {
		    uint32_t m = p/4;
		    mpz_set_ui(pol.get_mpz_t(), p);
		    mpz_powm_ui(temp.get_mpz_t(), temp.get_mpz_t(), m+1, pol.get_mpz_t());
		    return mpz_get_ui(temp.get_mpz_t());
	    } else { // p % 8 == 5
		    uint32_t m = p/8;
		    mpz_set_ui(pol.get_mpz_t(), p);
		    mpz_powm_ui(temp2.get_mpz_t(), temp.get_mpz_t(), 2*m+1, pol.get_mpz_t());
		    if (mpz_cmp_ui(temp2.get_mpz_t(), 1) == 0){
			    mpz_powm_ui(temp2.get_mpz_t(), temp.get_mpz_t(), m+1, pol.get_mpz_t());
			    return mpz_get_ui(temp2.get_mpz_t());
		    } else {
	    		mpz_mul_ui(temp.get_mpz_t(), temp.get_mpz_t(), 4);
		    	mpz_powm_ui(temp2.get_mpz_t(), temp.get_mpz_t(), m+1, pol.get_mpz_t());
    			if (mpz_divisible_ui_p(temp2.get_mpz_t(), 2)){
                    return mpz_get_ui(temp2.get_mpz_t())/2;
	    		} else {
                    return (mpz_get_ui(temp2.get_mpz_t())+p)/2;
			    }
		    }
	    }
    }

    uint32_t find_root_tonelli_shanks(uint32_t p, mpz_class &temp) {
        mpz_class pol;
        mpz_class temp2;
        uint32_t s = p-1;
        uint32_t e = 0;
        uint32_t a = mpz_get_ui(temp.get_mpz_t());
        while (s % 2 == 0){
            s /= 2;
            e ++;
        }

        mpz_set_ui(pol.get_mpz_t(), p);
        mpz_set_ui(temp2.get_mpz_t(), 2);
        mpz_powm_ui(temp.get_mpz_t(), temp2.get_mpz_t(), (p-1)/2, pol.get_mpz_t());
        while (mpz_cmp_ui(temp.get_mpz_t(), p-1) != 0){
            mpz_add_ui(temp2.get_mpz_t(), temp2.get_mpz_t(), 1);
            mpz_powm_ui(temp.get_mpz_t(), temp2.get_mpz_t(), (p-1)/2, pol.get_mpz_t());
        }
        
        mpz_class x;
        mpz_class b;
        mpz_class g;
        mpz_set_ui (temp.get_mpz_t(), a);
        mpz_powm_ui (x.get_mpz_t(), temp.get_mpz_t(), (s+1)/2, pol.get_mpz_t());
        mpz_powm_ui (b.get_mpz_t(), temp.get_mpz_t(), s, pol.get_mpz_t());
        mpz_powm_ui (g.get_mpz_t(), temp2.get_mpz_t(), s, pol.get_mpz_t());
        uint32_t r = e;
        
        while(true) {
            uint32_t m = 0;
            mpz_set (temp.get_mpz_t(), b.get_mpz_t());
            while (mpz_cmp_ui(temp.get_mpz_t(), 1) != 0) {
                mpz_mul(temp.get_mpz_t(), temp.get_mpz_t(), temp.get_mpz_t());
                mpz_mod(temp.get_mpz_t(), temp.get_mpz_t(), pol.get_mpz_t());
                m++;
            }
            if (m == 0) {
                return mpz_get_ui(x.get_mpz_t());
            } else {
                mpz_ui_pow_ui(temp.get_mpz_t(), 2, r-m-1);
                mpz_powm(temp.get_mpz_t(), g.get_mpz_t(), temp.get_mpz_t(), pol.get_mpz_t());
                mpz_mul(x.get_mpz_t(), x.get_mpz_t(), temp.get_mpz_t());
                mpz_mod(x.get_mpz_t(), x.get_mpz_t(), pol.get_mpz_t());
                mpz_mul(temp.get_mpz_t(), temp.get_mpz_t(), temp.get_mpz_t());
                mpz_mod(temp.get_mpz_t(), temp.get_mpz_t(), pol.get_mpz_t());
                mpz_mul(b.get_mpz_t(), b.get_mpz_t(), temp.get_mpz_t());
                mpz_mod(b.get_mpz_t(), b.get_mpz_t(), pol.get_mpz_t());
                mpz_set(g.get_mpz_t(), temp.get_mpz_t());
                r = m;
            }
        }
    }

    uint32_t find_root (mpz_class &k, uint32_t p) { // finds modular square root of k (mod p)
        mpz_class temp;
	    mpz_mod_ui(temp.get_mpz_t(), k.get_mpz_t(), p);
	    uint32_t a = mpz_get_ui(temp.get_mpz_t());
        if (p == 2) {
            return a;
        }
        else if (p % 4 == 3 || p % 8 == 5) {
            return find_root_pocklington(p, temp);
        }
        else {
            return find_root_tonelli_shanks(p, temp);
        }
    }

    void initialize(int p1, int p2, double f) {
        fb_bound = static_cast<uint32_t>(PARAMS_ARR[p1].FBB * f + PARAMS_ARR[p2].FBB * (1 - f));
        lp_bound = fb_bound * static_cast<uint32_t>(PARAMS_ARR[p1].LPB * f + PARAMS_ARR[p2].LPB * (1 - f));
        M = static_cast<uint32_t>(PARAMS_ARR[p1].M * f + PARAMS_ARR[p2].M * (1 - f));
        T = static_cast<float>(PARAMS_ARR[p1].T * f + PARAMS_ARR[p2].T * (1 - f));

        int max_idx = -1;
        double max_score = get_multiplier_score(1);
        for (int i = 0; i < 18; i++) {
            double score = get_multiplier_score(SMALL_PRIMES[i]);
            if (score > max_score) {
                max_score = score;
                max_idx = i;
            }
        }

        multiplier = (max_idx == -1) ? 1 : SMALL_PRIMES[max_idx];
        if (max_idx != -1) {
            mpz_mul_ui(nc.get_mpz_t(), nc.get_mpz_t(), SMALL_PRIMES[max_idx]);
        }
        
        std::vector<char> vals = era_sieve(fb_bound);
        fb = extract_primes(fb_bound, multiplier, vals);
        roots.resize(fb.size());
        fb_logs.resize(fb.size());
        for (size_t i = 0; i < fb.size(); i++) {
            if (fb[i] == multiplier) {
                roots[i] = 0;
            }
            else {
                roots[i] = find_root(nc, fb[i]);
                if (roots[i] > fb[i]/2) {
                    roots[i] = fb[i] - roots[i];
                }
            }
            fb_logs[i] = fast_log(fb[i]);
        }

        // TODO allocate space for matrix relations
        // TODO initialize partials
    }

    Sieve(mpz_class &nc) : nc(nc)
    {
        int bs = mpz_sizeinbase(nc.get_mpz_t(), 2);

        if (bs <= PARAMS_ARR[0].bits) {
            initialize(0, 0, 0);
        }
        else if (bs >= PARAMS_ARR[9].bits) {
            initialize(9, 9, 0);
        }
        else {
            auto it = std::find_if(PARAMS_ARR.begin(), PARAMS_ARR.end(), [&bs](auto &a) { return a.bits >= bs; });
            size_t i = std::distance(PARAMS_ARR.begin(), it);
            initialize(i, i-1, (bs - PARAMS_ARR[i-1].bits) / (PARAMS_ARR[i].bits - PARAMS_ARR[i-1].bits));
        }
    }
};
