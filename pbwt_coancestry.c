#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pbwtmaster.h"

int pbwt_coancestry(const cmd_t *c)
{
    int v = 0;
    size_t i = 0;
    size_t j = 0;
    pbwt_t *b = NULL;

    if (c == NULL)
    {
        return -1;
    }

    /* Read in the pbwt file from disk */
    b = pbwt_read(c->instub);
    if (b == NULL)
    {
        fprintf(stderr, "pbwtmaster [ERROR]: cannot read data from %s\n", c->instub);
        return -1;
    }

    /* Uncompress the haplotype data */
    v = pbwt_uncompress(b);
    if (v < 0)
    {
        fputs("pbwtmaster [ERROR]: error uncompressing haplotype data\n", stderr);
        return -1;
    }

    /* Construct adjacency list */
    if (c->adjlist)
    {
        if (c->set_match)
        {
            if (c->print_sites)
            {
                v = pbwt_set_match(b, c->minlen, report_adjlist_with_sites);
            }
            else
            {
                v = pbwt_set_match(b, c->minlen, report_adjlist);
            }
        }
        else
        {
            if (c->print_sites)
            {
                v = pbwt_all_match(b, c->minlen, report_adjlist_with_sites);
            }
            else
            {
                v = pbwt_all_match(b, c->minlen, report_adjlist);
            }
        }
    }
    else if (c->count_only)
    {
        if (c->set_match)
        {
            v = pbwt_set_match(b, c->minlen, add_nmatch);
        }
        else
        {
            v = pbwt_all_match(b, c->minlen, add_nmatch);
        }
        if (c->out_diploid)
        {
            for (i = 0; i < b->nsam/2; ++i)
            {
                for (j = 0; j < b->nsam/2 - 1; ++j)
                {
                    printf("%zu\t", b->nmatrix[2*i][2*j] + b->nmatrix[2*i+1][2*j+1]);
                }
                printf("%zu\n", b->nmatrix[2*i][2*j+1] + b->nmatrix[2*i+1][2*j]);
            }
        }
        else
        {
            /* Print coancestry matrix to STDOUT */
            for (i = 0; i < b->nsam; ++i)
            {
                for (j = 0; j < b->nsam - 1; ++j)
                {
                    printf("%zu\t", b->nmatrix[i][j]);
                }
                printf("%zu\n", b->nmatrix[i][j]);
            }
        }
    }
    else
    {
        /* Find matches */
        if (c->set_match)
        {
            v = pbwt_set_match(b, c->minlen, add_coancestry);
        }
        else
        {
            v = pbwt_all_match(b, c->minlen, add_coancestry);
        }
        if (v < 0)
        {
            fputs("pbwtmaster [ERROR]: error retrieving matches\n", stderr);
            return -1;
        }

        if (c->out_diploid)
        {
            for (i = 0; i < b->nsam/2; ++i)
            {
                for (j = 0; j < b->nsam/2 - 1; ++j)
                {
                    printf("%1.4lf\t", b->cmatrix[2*i][2*j] + b->cmatrix[2*i+1][2*j+1]);
                }
                printf("%1.4lf\n", b->cmatrix[2*i][2*j+1] + b->cmatrix[2*i+1][2*j]);
            }
        }
        else
        {
            /* Print coancestry matrix to STDOUT */
            for (i = 0; i < b->nsam; ++i)
            {
                for (j = 0; j < b->nsam - 1; ++j)
                {
                    printf("%1.4lf\t", b->cmatrix[i][j]);
                }
                printf("%1.4lf\n", b->cmatrix[i][j]);
            }
        }
    }

    /* Clean up allocated memory */
    pbwt_destroy(b);

    return 0;
}
