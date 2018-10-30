# libplink

An optimized and light-weight C library for reading and writing plink format data

## How to build

First clone the git repository

```
git clone https://github.ancestry.com/dgarrigan/libplink.git
```

Then change directory and use the Makefile provided

```
cd libplink; make all
```

This will create a static library file `libplink.a`, which can be linked when you build your executable.
Included in the distribution is an example program that uses the libplink API. To build this program, type

```
gcc -Wall -O2 -L. -o test-api test-api.c -lplink
```

to create the executable `test-api`.


## Using the API

The API currently provides three basic functions to the user

```
plink_t * read_plink (const char *instub, const int has_reg, const int is_phased)
```

which reads all data from plink set into memory as a plink_t object (see below). The function
needs to know if there is an associated .reg file with the plink data set and whether there the
data are unphased (.bed) or phased (.hap).

```
char * hap2str (plink_t *p, const uint64_t i, const int parent)
```

which retrieves one haplotype from sample *i* and returns it as a null-terminated string. The parent parameter is 0 for the first haplotype or 1 to retreive the second haplotype.

```
char * query_reg (reg_t *r, const char *iid)
```

This function queries the data loaded from a .reg file. The query term is the sample id and the function returns a string of the region name associated with that sample in the .reg file.

## Contact

For questions or bug reports, contact [daniel.garrigan@ancestry.com](mailto:daniel.garrigan@ancestry.com)