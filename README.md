disambiguate
============

Disambiguation algorithm for reads aligned to human and mouse genomes using Tophat or BWA mem.
Both a Python and C++ implementation will be offered. The Python implementation has a dependency
on the Pysam module. The C++ implementation will depend on the availability of zlib and the
Bamtools C++ API.

For usage help, run disambiguate.py as-is.

To compile the C++ program, use the following syntax:
c++ -I /path/to/bamtools_c_api/include/ -L /path/to/bamtools_c_api/lib/ -o disambiguate dismain.cpp -lz -lbamtools
