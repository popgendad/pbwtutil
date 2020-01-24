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
            v = pbwt_set_match(b, c->minlen, report_adjlist);
        }
        else
        {
            v = pbwt_all_match(b, c->minlen, report_adjlist);
        }
        /* TODO: Implement adjlist diploid
        if (c->out_diploid)
        {
            adjlist_t *h = diploidize(g);
            print_adjlist(h);
        }
        else
        {
            print_adjlist(g);
        }*/
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
                    printf("%1.4lf\t", b->cmatrix[2*i][2*j+1] + b->cmatrix[2*i+1][2*j]);
                }
                printf("%1.4lf\n", b->cmatrix[2*i][2*j] + b->cmatrix[2*i+1][2*j+1]);
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
