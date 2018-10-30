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

## Data structures

The main `plink_t` structure is defined this way

```
typedef struct _plink_t
{
    size_t nsam;    /* Number of diploid samples */
    size_t nsnp;    /* Number of single nucleotide polymorphism loci */
    bed_t *bed;     /* Pointer to array of bed data structures */
    bim_t *bim;     /* Pointer to array of bim data structures */
    fam_t *fam;     /* Pointer to array of fam data structures */
    reg_t *reg;     /* Pointer to array of reg data structures */
} plink_t;
```

Likewise, there are data structures defined for each component. The `bim_t` structure is

```
typedef struct _bim_t
{
    int ch;           /* plink-format integer chromosome identifier */
    char *rsid;       /* marker ID, generally matches '^rs[0-9]+$' */
    double cM;        /* recombination distance from previous marker (or beginning of chromosome) */
    uint64_t bp;      /* coordinate location in base pairs */
    char a0;          /* First allele in genotype */
    char a1;          /* Second allele in genotype */
    khash_t(integer) *index;  /* Hash keyed on rsid */
} bim_t;
```

The `fam_t` data structure

```
typedef struct _fam_t
{
    char *fid;           /* Family ID */
    char *iid;           /* Individual ID */
    char *pid;           /* Paternal ID */
    char *mid;           /* Maternal ID */
    char *sex;           /* Sex */
    char *phe;           /* Encoded phenotype */
    khash_t(integer) *index;  /* Hash keyed on iid */
} fam_t;
```

The reg data structure `reg_t`

```
typedef struct _reg_t
{
    char *fid;           /* Family ID */
    char *iid;           /* Individual ID */
    char *pid;           /* Paternal ID */
    char *mid;           /* Maternal ID */
    char *sex;           /* Sex */
    char *phe;           /* Encoded phenotype */
    char *pop;           /* Population name */
    char *reg;           /* Region name */
    khash_t(integer) *index;  /* Hash keyed on iid */
} reg_t;
```

And finally the bed/hap files, which hold all of the genotype data, type `bed_t`

```
typedef struct _bed_t
{
    int header1;
    int header2;
    int phased;
    int orientation;
    uint64_t record_size;       /* number of bytes per record (major-order singleton array) */
    uint64_t size;              /* number of bytes of data */
    unsigned char *data;
} bed_t;
```

## Contact

For questions or bug reports, contact [daniel.garrigan@ancestry.com](mailto:daniel.garrigan@ancestry.com)