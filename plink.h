#ifndef PLINK_H
#define PLINK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "khash.h"

/************************************************
 * Definitions
 ************************************************/

#define VERSION 0.0.3-alpha

/* Hash types */
KHASH_MAP_INIT_STR(string, char *)
KHASH_MAP_INIT_STR(integer, uint64_t)

/* BED orientation */
#define INDIVIDUAL_MAJOR_ORDER 0
#define SNP_MAJOR_ORDER 1

/* Two header bytes for bed/hap data */
#define BED_MAGIC1 108
#define BED_MAGIC2 27
#define HAP_MAGIC1 108
#define HAP_MAGIC2 211
#define HAP_00 0
#define HAP_01 1
#define HAP_10 2
#define HAP_11 3

/* Allele encodings */
#define BED_HOMOZYGOUS_0 0
#define BED_MISSING 1
#define BED_HETEROZYGOUS 2
#define BED_HOMOZYGOUS_1 3

/************************************************
 * Globally scoped variables
 ************************************************/


/************************************************
 * Data structures
 ************************************************/

/* SNP entry from BIM file */
typedef struct _bim_t
{
    int ch;           /* chromosome number--kept as int for sorting purposes */
    char *rsid;       /* marker ID--generally matches '^rs[0-9]+$' */
    double cM;        /* recombination distance from previous marker (or beginning of chromosome) */
    uint64_t bp;      /* coordinate location in base pairs */
    char a0;          /* e.g., 'A','G'; represented by 0 and 1 in a .bed (or .hap) file, respectively (typically allele0 would be the major allele, but this may not be guaranteed).  Set to ALLELE_MISSING if the allele value is as of yet unknown. */
    char a1;
    khash_t(integer) *index;
} bim_t;

/* Individual entry from FAM file */
typedef struct _fam_t
{
    char *fid;
    char *iid;
    char *pid;
    char *mid;
    char *sex;
    char *phe;
    khash_t(integer) *index;
} fam_t;

/* Individual entry from REG file */
typedef struct _reg_t
{
    char *fid;
    char *iid;
    char *pid;
    char *mid;
    char *sex;
    char *phe;
    char *pop;
    char *reg;
    khash_t(integer) *index;
} reg_t;

/* Binary genotype data and associated functionality */
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

/* Main plink_t class definition */
typedef struct _plink_t
{
    size_t nsam;
    size_t nsnp;
    bed_t *bed;
    bim_t *bim;
    fam_t *fam;
    reg_t *reg;
} plink_t;

/************************************************
 * Function prototypes
 ************************************************/

/* Read input bed/hap file.  Use given unsigned char array or allocate a new one if argument is NULL */
extern bed_t *read_bed (const char *, uint64_t n_indiv, uint64_t n_snps, unsigned char *data);

/* Write bed/hap to file, return number of bytes written */
extern uint64_t write_bed (FILE *, const bed_t *);

/* Read input bim file */
extern bim_t *read_bim (const char *, size_t *);

/* Index a bim dataset */
extern khash_t(integer) *index_bim (const bim_t *, const size_t);

/* Write marker information to bim file */
extern int write_bim (const char *, const bim_t *, const size_t);

/* Read input fam file */
extern fam_t *read_fam (const char *, size_t *);

/* Index a fam dataset */
extern khash_t(integer) *index_fam (const fam_t *, const size_t);

/* Write sample information to fam file */
extern int write_fam (const char *, const fam_t *, const size_t);

/* Read input reg file */
extern reg_t *read_reg (const char *, size_t *);

/* Index a reg data set */
extern khash_t(integer) *index_reg (const reg_t *, const size_t);

/* Write region information to reg file */
extern int write_reg (const char *, const reg_t *, const size_t);


/************************************************
 * API function prototypes
 ************************************************/

/* Read all data from plink set into memory */
extern plink_t *read_plink (const char *, const int, const int);

/* Get haplotype string */
extern char *hap2str (plink_t *, const uint64_t, const int);

/* Query the reg data with a sample id and get a region name back */
extern char *query_reg (reg_t *, const char *);

/* Get array of unsigned long integers representing binary encoding of SNPs */
extern uint64_t *hap2ulong (plink_t *, const uint64_t, const int);


/************************************************
 * Inline functions
 ************************************************/

/* Compute bed data binary record size based on number of elements in each record
 * (e.g., if using SNP-major order, this is the number of inidiviuals) */
static inline uint64_t
bed_record_size (int n_minor)
{
    return n_minor / 4 + ((n_minor % 4) ? 1 : 0);
}

/* retreive 2-bit encoding for the given major (e.g., marker index if we're doing SNP-major order) and minor (e.g., individual index) matrix entry */
static inline unsigned char
get_encoding (const unsigned char *data, uint64_t record_size, uint64_t major, uint64_t minor)
{
    return (data[major * record_size + minor / 4] >> 2 * (minor % 4)) & 3;
}

/* Get unphased genotype */
static inline unsigned char
genotype (bed_t *bed, uint64_t i, uint64_t m)
{
    return get_encoding(bed->data, bed->record_size, bed->orientation ? m : i, bed->orientation ? i : m);
}

/* Get phased genotype */
static inline int
plink_haplotype (bed_t *bed, uint64_t i, uint64_t m, int parent)
{
    return (genotype(bed, i, m) & (1 << parent)) ? 1 : 0;
}

#endif