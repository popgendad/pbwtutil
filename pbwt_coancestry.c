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
        fputs("pbwtmaster [ERROR]: error uncompressing haplotype data", stderr);
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
/*        if (c->out_diploid)
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
        if (c->out_diploid)
        {
/*            size_t new_nsam = b->nsam / 2;


            cmatrix = (double **)malloc(new_nsam * sizeof(double *));
            if (cmatrix == NULL)
            {
                fputs("pbwtmaster [ERROR]: memory allocation error", stderr);
                return -1;
            }
            for (i = 0; i < new_nsam; ++i)
            {
                cmatrix[i] = (double *)malloc(new_nsam * sizeof(double));
                if (cmatrix[i] == NULL)
                {
                    fputs("pbwtmaster [ERROR]: memory allocation error", stderr);
                    return -1;
                }
                memset(cmatrix[i], 0, new_nsam * sizeof(double));
            }


            if (c->set_match)
            {
                v = pbwt_set_match(b, c->minlen);
            }
            else
            {
                v = pbwt_all_match(b, c->minlen, report_match);
            }
            if (v < 0)
            {
                fputs("pbwtmaster [ERROR]: error retrieving matches", stderr);
                return -1;
            }


            match_coasearch(b, b->match, cmatrix, 0, b->nsite, c->out_diploid);


            for (i = 0; i < new_nsam; ++i)
            {
                for (j = 0; j < new_nsam - 1; ++j)
                {
                    printf("%1.4lf\t", cmatrix[i][j]);
                }
                printf("%1.4lf\n", cmatrix[i][j]);
            }

            for (i = 0; i < new_nsam; ++i)
            {
                free(cmatrix[i]);
            }
            free(cmatrix);*/
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
                fputs("pbwtmaster [ERROR]: error retrieving matches", stderr);
                return -1;
            }

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
