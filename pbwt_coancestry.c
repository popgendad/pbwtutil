#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ancmatch.h"

int pbwt_coancestry(cmd_t *c)
{
    int v = 0;
    size_t i = 0;
    size_t j = 0;
    size_t nsites = 0;
    double **cmatrix = NULL;
    pbwt_t *b = NULL;

    /* Read in the pbwt file from disk */
    b = pbwt_read(c->instub);
    if (b == NULL)
    {
        fprintf(stderr, "ancmatch [ERROR]: cannot read data from %s\n", c->instub);
        return 1;
    }

	/* Unless specified otherwise, indicate penultimate site is end site */
    nsites = b->nsite;

    /* Uncompress the haplotype data */
    pbwt_uncompress(b);

    /* Initialize coancestry matrix */
    cmatrix = (double **)malloc(b->nsam * sizeof(double *));
    for (i = 0; i < b->nsam; ++i)
    {
        cmatrix[i] = (double *)malloc(b->nsam * sizeof(double));
        for (j = 0; j < b->nsam; ++j)
        {
            cmatrix[i][j] = 0.0;
        }
    }

    /* Find matches to fill coancestry matrix */
    v = pbwt_set_match(b, c->minlen);
    match_search(b, b->match, cmatrix, 0, nsites);
    for (i = 0; i < b->nsam; ++i)
    {
        for (j = 0; j < b->nsam; ++j)
        {
            if (j < b->nsam - 1)
            {
                printf("%1.4lf\t", cmatrix[i][j]);
            }
            else
            {
                printf("%1.4lf\n", cmatrix[i][j]);
            }
        }
    }

    /* Clean up allocated memory */
    pbwt_destroy(b);

    return v;
}
