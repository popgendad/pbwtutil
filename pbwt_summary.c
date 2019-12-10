#include <stdio.h>
#include <stdlib.h>
#include "pbwtmaster.h"


int pbwt_summary(cmd_t *c)
{
    int v = 0;
    size_t orig_size = 0;
    pbwt_t *b = NULL;

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


    /* Print summary report */
    printf("Number of samples:\t%zu\n", b->nsam);
    printf("Number of sites:\t%zu\n", b->nsite);
    printf("Total recombination distance:\t%1.5lf\n", b->cm[b->nsite-1] - b->cm[0]);
    printf("Size of compressed data:\t%zu\n", orig_size);
    printf("Size of uncompressed data:\t%zu\n", b->datasize);

    /* Clean up allocated memory */
    pbwt_destroy(b);
    free(c->query);
    free(c);

    return v;
}