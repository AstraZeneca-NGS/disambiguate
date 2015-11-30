disambiguate
============

Disambiguation algorithm for reads aligned to two species (e.g. human and mouse genomes) from 
Tophat, Hisat2, STAR or BWA mem. Both a Python and C++ implementation are offered. The Python 
implementation has a dependency on the Pysam module. The C++ implementation depends on the 
availability of zlib and the Bamtools C++ API. For STAR alignments it is highly recommended
to include the NM tag in the output when performing alignment (in fact this is a requirement
for the C++ version).

Differences between the Python and C++ versions:
1. The Python version can do name sorting of the reads (a necessary step) internally
but for the C++ version the input BAM files _must_ be name sorted (internal name sorting not 
supported).
2. The flag -s (samplename prefix) must be provided as an input parameter to the C++ binary

For usage help, run disambiguate.py as-is.

To compile the C++ program, use the following syntax in the same folder where the code is:
c++ -I /path/to/bamtools_c_api/include/ -I./ -L /path/to/bamtools_c_api/lib/ -o disambiguate dismain.cpp -lz -lbamtools
