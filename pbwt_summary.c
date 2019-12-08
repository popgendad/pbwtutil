#include <stdio.h>
#include <stdlib.h>
#include "pbwtmaster.h"


int pbwt_summary(cmd_t *c)
{
    int v = 0;
    size_t orig_size = 0;
    pbwt_t *b = NULL;

    b = pbwt_read(c->instub);
    orig_size = b->datasize;
    v = pbwt_uncompress(b);
    printf("Number of samples:\t%zu\n", b->nsam);
    printf("Number of sites:\t%zu\n", b->nsite);
    printf("Size of compressed data:\t%zu\n", orig_size);
    printf("Size of uncompressed data:\t%zu\n", b->datasize);
    pbwt_destroy(b);

    return v;
}