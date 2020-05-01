Mash is normally distributed as a dependency-free binary for Linux or OSX (see
https://github.com/marbl/Mash/releases). This source distribution is intended
for other operating systems or for development. Mash requires c++14 to build,
which is available in and GCC >= 5 and XCode >= 6.

See http://mash.readthedocs.org for more information.

# TODO
- [x] fix missing hash value
- [x] add bad check
- [x] optimize screen operation
- [ ] support fast fastq parser to screen
- [ ] support gziped files
- [ ] fix addMinHash for sketch whole file when checking bad characters
- [ ] free mermory intermediately for memory limit
- [x] bitwise reverseComplement doesn't fit unclean databases.
- [ ] the SIMD MurmurHash3 can't work normally when kmerSize is greater than 64.
- [ ] the MurmurHash3_x86_128 version isn't vectorized.
