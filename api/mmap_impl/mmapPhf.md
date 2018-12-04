# Introduction #

This is a simple implementation of the CHD perfect hash algorithm. CHD can
generate perfect hash functions for very large key sets--on the order of
millions of keys--in a very short time. On my circa 2012 desktop and using
the default parameters (hash load factor of 80% and average displacement map
bucket load of 4.0 keys) this implementation can generate a hash function
for 1,000 keys in less than 1/100th of a second, and 1,000,000 keys in less
than a second.

This code was taken PHF package by William Ahern,
https://github.com/wahern/phf for use in HAL.
It has been modified from the original source.  Copied at commit
2b781079946401ce160eb577b4534e94b763242a.


## C++ API ##

### PHF::uniq<T>(T k[], size_t n); ###

Similar to the shell command `sort | uniq`. Sorts, deduplicates, and shifts
down the keys in the array k. Returns the number of unique keys, which will
have been moved to the beginning of the array. If necessary do this before
calling PHF::init, as PHF::init does not tolerate duplicate keys.

### int PHF::init<T, nodiv>(struct phf *f, const T k[], size_t n, size_t l, size_t a, phf_seed_t s);

Generate a perfect hash function for the n keys in array k and store the
results in f. Returns a system error number on failure, or 0 on success. f
is unmodified on failure.

* keys: array of keys in order from 1..#keys. They should be all
  numbers or all strings.

* lambda: number of keys per bucket when generating the g() function mapping.

* alpha: output hash space loading factor as percentage from
  1..100. 100% generates a *minimal* perfect hash function. But note that
  the implementation does *not* implement the necessary optimizations to
  ensure timely generation of minimal perfect hash functions. Normally you
  want a loading factor of 80% to 90% for large key sets.

* seed: random integer seed.

* nodiv: if true rounds r and m to powers of 2, and performs modular
  reduction using bitwise AND. Otherwise, r and m are rounded up to the
  nearest primes and modulo division used when indexing tables. Note that
  the rounding occurs after calculation of the intermediate and output hash
  table loading.

  This is more important when building small hash tables with the C
  interface. The optimization is substantial when the compiler can inline
  the code, but isn't substantial from Lua.

### void PHF::destroy(struct phf *);

Deallocates internal tables, but not the struct object itself.

### void PHF::compact<T, nodiv>(struct phf *);

By default the displacement map is an array of uint32_t integers. This
function will select the smallest type necessary to hold the largest
displacement value and update the internal state accordingly. For a loading
factor of 80% (0.8) in the output hash space, and displacement map loading
factor of 4 (400%), the smallest primitive type will often be uint8_t.

### phf_hash_t PHF::hash<T>(struct phf *f, T k);

Returns an integer hash value, h, where 0 <= h < f->m. h will be unique for
each unique key provided when generating the function. f->m will be larger
than the number of unique keys and is based on the specified loading factor
(alpha), rounded up to the nearest prime or nearest power of 2, depending on
the mode of modular reduction selected. For example, for a loading factor of
80% m will be 127: 100 is 80% of 125, and 127 is the closest prime greater
than or equal to 125. With the nodiv option, m would be 128: 100 is 80% of
125, and 128 is the closest power of 2 greater than or equal to 125.

