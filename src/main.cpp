#include <boost/program_options.hpp>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <gmpxx.h>
#include <gmp.h>
#include <mutex>
#include "sieve.h"

using namespace boost;
namespace po = boost::program_options;

void print_factors(const std::vector<std::string> &factors) {
    if (factors.size() == 2) { // expected
        std::cout << "p: " << factors[0] << std::endl;
        std::cout << "q: " << factors[1] << std::endl;
    }
    else { // unexpected: display all factors anyway
        for (auto factor : factors) {
            std::cout << factor << std::endl;
        }
    }
}

std::vector<std::string> trial_division(mpz_class &nc, int bound)
{
    std::vector<std::string> factors;

    while (mpz_divisible_ui_p (nc.get_mpz_t(), 2)){
		mpz_divexact_ui (nc.get_mpz_t(), nc.get_mpz_t(), 2);
		factors.push_back(std::to_string(2));
	}

    long t = 3;
	while (t < bound){
		while (mpz_divisible_ui_p (nc.get_mpz_t(), t)){
			mpz_divexact_ui (nc.get_mpz_t(), nc.get_mpz_t(), t);
            factors.push_back(std::to_string(t));
		}
		t += 2;
	}

    return factors;
}

std::vector<std::string> pollards_rho(mpz_class &nc, int steps)
{
    std::vector<std::string> factors;

    mpz_class xc;
    mpz_class yc;
    mpz_class gc;
    mpz_class tempc;

    int ct = 0;
    while (ct < steps) {
        mpz_set_ui(gc.get_mpz_t(), 1);
        for (int i = 0; i < 128; i++) {
            mpz_mul(tempc.get_mpz_t(), xc.get_mpz_t(), xc.get_mpz_t());
            mpz_add_ui(tempc.get_mpz_t(), tempc.get_mpz_t(), 1);
            mpz_mod(xc.get_mpz_t(), tempc.get_mpz_t(), nc.get_mpz_t());
            mpz_mul(tempc.get_mpz_t(), yc.get_mpz_t(), yc.get_mpz_t());
            mpz_add_ui(tempc.get_mpz_t(), tempc.get_mpz_t(), 1);
            mpz_mul(yc.get_mpz_t(), tempc.get_mpz_t(), tempc.get_mpz_t());
            mpz_add_ui(tempc.get_mpz_t(), yc.get_mpz_t(), 1);
            mpz_mod(yc.get_mpz_t(), tempc.get_mpz_t(), nc.get_mpz_t());
            mpz_sub(tempc.get_mpz_t(), xc.get_mpz_t(), yc.get_mpz_t());
            mpz_mul(gc.get_mpz_t(), gc.get_mpz_t(), tempc.get_mpz_t());
            mpz_mod(gc.get_mpz_t(), gc.get_mpz_t(), nc.get_mpz_t());
        }

        mpz_gcd(gc.get_mpz_t(), gc.get_mpz_t(), nc.get_mpz_t());
        if (mpz_cmp_ui(gc.get_mpz_t(), 1) != 0) {
            if (mpz_probab_prime_p(gc.get_mpz_t(), 12)) {
                factors.push_back(gc.get_str());
                mpz_divexact(nc.get_mpz_t(), nc.get_mpz_t(), gc.get_mpz_t());
            }
            else { // composite
                factors.push_back(gc.get_str());
                mpz_divexact(nc.get_mpz_t(), nc.get_mpz_t(), gc.get_mpz_t());
            }
            if (mpz_cmp_ui(nc.get_mpz_t(), 1) == 0) {
                break;
            }
            if (mpz_probab_prime_p(nc.get_mpz_t(), 12)) {
                factors.push_back(nc.get_str());
                break;
            }
        }

        ct += steps;
    }

    return factors;
}

void *run_sieve_thread(void *args) {

    // TODO call polygroup_init
    // TODO call generate_polygroup
    // TODO call poly_init
    // TODO call generate_poly
    // TODO call sieve_poly
    // TODO call add_polygroup_relations
    // TODO update mutex data

    return nullptr;
}

void multithreaded_factor(Sieve &sv, int threads) {

    // TODO call gpool_init
    // TODO call advance_gpool
    // TODO call pthread_create (run_sieve_thread)
    // TODO call pthread_join
    // TODO call build_matrix
    // TODO call solve_matrix
}

std::vector<std::string> quadratic_sieve(mpz_class &nc, int steps)
{
    std::vector<std::string> factors;

    Sieve sieve(nc);

    // TODO call multithreaded factor

    return factors;
}

int main(int argc, char *argv[])
{
    try {
        int threads;
        std::string number; 
        po::options_description options("Allowed options");
        options.add_options()
            ("help", "produce help message")
            ("threads,t", po::value<int>(&threads)->default_value(1), "number of threads")
            ("number,n", po::value<std::string>(&number), "number to factor")
        ;

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(options).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << "Usage: prime_sieves [options]" << std::endl;
            std::cout << options;
            return 0;
        }

        std::cout << "-------------\nprime_sieves\n-------------" << std::endl;
        std::cout << "starting prime sieve for: " << number << std::endl;
        std::cout << "using " << threads << " thread(s)\n" << std::endl;

        // convert number to mpz_t
        mpz_class nc;
        nc = number.c_str();

        std::cout << "searching for small factors of ``n`` using trial division..." << std::endl;
        std::vector<std::string> factors = trial_division(nc, 32768);
        if (!factors.empty()) {
            print_factors(factors);
            std::cout << "finished (solved using trial division)\n\n";
            return 0;
        }

        std::cout << "searching for small factors of ``n`` using pollards rho..." << std::endl;
        factors = pollards_rho(nc, 65536);
        if (!factors.empty()) {
            print_factors(factors);
            std::cout << "finished (solved using pollards rho)\n\n";
            return 0;
        }

        std::cout << "searching for whether ``n`` is a prime number..." << std::endl;
        if (mpz_cmp_ui(nc.get_mpz_t(), 1) == 0 || mpz_probab_prime_p(nc.get_mpz_t(), 15)) {
            std::cout << "1" << std::endl;
            std::cout << nc.get_str() << std::endl;
            std::cout << "finished (n is prime)\n\n";
            return 0;
        }

        std::cout << "searching for large factors of ``n`` using quadratic sieve..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        factors = quadratic_sieve(nc, threads);
        auto end = std::chrono::high_resolution_clock::now(); 
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); 
        std::cout << "finished quadratic sieve in " << duration.count() << "ms\n\n";
    }
    catch(std::exception& e)
    {
        std::cout << e.what() << "\n";
        return 1;
    }

    return 0;
}