# prime_sieves

Implementation of Quadratic Sieve.

## Dependencies

Libraries used:

* gmp
* pthread
* boost

<!-- install using sudo apt-get install libgmp-dev -->

## Instructions

Use following commands to build and run program on linux:

```bash
mkdir build
cmake .. && make

# or whatever number you want
./main 15347 -threads 4
```

## Sample Output

n = 15437

```bash
-------------
prime_sieves
-------------
starting prime sieve for: 15437
using 1 thread(s)

searching for small factors of ``n`` using trial division...
p: 43
q: 359
finished (solved using trial division)
```

n = 1750793698402698365078886617619033336333347619191

```bash
# TODO
```