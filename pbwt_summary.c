#include <stdio.h>
#include <stdlib.h>
#include "pbwtmaster.h"


int pbwt_summary(const cmd_t *c)
{
    int v = 0;
    size_t orig_size = 0;
    khash_t(integer) *regcount = NULL;
    khint_t it = 0;
    pbwt_t *b = NULL;

    if (c == NULL)
    {
        return -1;
    }

    /* Read PBWT file into memory */
    b = pbwt_read(c->instub);
    if (b == NULL)
    {
        fprintf(stderr, "pbwtmaster [ERROR]: cannot read data from %s\n", c->instub);
        return -1;
    }

    /* Get size of compressed haplotype data */
    orig_size = b->datasize;

    /* Uncompress haplotype data */
    v = pbwt_uncompress(b);
    if (v < 0)
    {
        fputs("pbwtmaster [ERROR]: error uncompressing haplotype data", stderr);
        return -1;
    }

    regcount = pbwt_get_regcount(b);

    /* Print summary report */
    printf("Number of samples:\t%zu\n", b->nsam);
    printf("Number of sites:\t%zu\n", b->nsite);
    printf("Total recombination distance:\t%1.5lf\n", b->cm[b->nsite-1] - b->cm[0]);
    printf("Number of regions:\t%d\n", kh_size(regcount));
    printf("Size of compressed data:\t%zu\n", orig_size);
    printf("Size of uncompressed data:\t%zu\n", b->datasize);

    /* If user specifies regcount option */
    if (c->reg_count)
    {
        putchar('\n');
        printf("Region\tCount\n");
        for (it = kh_begin(regcount); it != kh_end(regcount); ++it)
        {
            if (kh_exist(regcount, it))
            {
                printf("%s\t%lu\n", kh_key(regcount, it), kh_value(regcount, it));
            }
        }
    }

    /* Clean up allocated memory */
    pbwt_destroy(b);

    return v;
}