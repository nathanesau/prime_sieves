# prime_sieves

Implementation of Prime Sieve Algorithms. I aim to implement the Quadratic Sieve algorithm in under 1000 lines of C++ code.

Libraries used:

* gmp.h

<!-- install using sudo apt-get install libgmp-dev -->

* pthread.h

## Test Cases

Test cases are in ``test_cases`` folder. 

* Test cases were created using primes in https://primes.utm.edu/curios/home.php.

## References

1. [nsieve](https://github.com/albinoBandicoot/nsieve/blob/master/README.technical) README provides some excellent details on prime sieve algorithms.

<!--
to compile, clone nsieve and on ubuntu run "sudo apt-get install libgmp3-dev"
then run make

results (these were maxing out my CPU usage at 100%):

    time ./nsieve 15347 # 0.012s
    time ./nsieve 1750793698402698365078886617619033336333347619191 -theads 12 # 0.116s
    time ./nsieve 1111114111144441147411444411114110934332856328033857561580338632856333351 -threads 12 # 16.612s
    time ./nsieve 798559766371214180866081220559317259672155849972506231361450196240562772751122717659661579 -threads 12 # 73m 53.662s

additional info for 90 digit case:

    Have 12822 of 15538 relations (6452 full + 6370 combined from 76909 partial); sieved 32379904 polynomials from 31621 groFatal error: Out of polynomials!!!

    num relations: 15538
    polynomials sieved: ~ 32,379,904
    msieve time: 26m 01s

    note: https://github.com/albinoBandicoot/nsieve/blob/master/src/rho.h is too slow.
    note: https://github.com/tuomas56/qsieve/blob/master/src/main.rs is too slow.

reference: nsieve is about 2000 lines of code.
-->