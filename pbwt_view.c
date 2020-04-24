#include <stdio.h>
#include <stdlib.h>
#include "pbwtmaster.h"

int pbwt_print(const pbwt_t *, const int);
int print_sites(const pbwt_t *);

int pbwt_view(const cmd_t *c)
{
    int v = 0;
    pbwt_t *b = NULL;

    if (c == NULL)
    {
        return -1;
    }

    /* Read PBWT file data into memory */
    b = pbwt_read(c->instub);
    if (b == NULL)
    {
        fprintf(stderr, "pbwtmaster [ERROR]: cannot read data from %s\n", c->instub);
        return -1;
    }

    /* Uncompress haplotype data */
    v = pbwt_uncompress(b);
    if (v < 0)
    {
        fputs("pbwtmaster [ERROR]: error uncompressing haplotype data\n", stderr);
        return -1;
    }

    /*ppa = pbwt_build(b); */

    /* Print the PBWT data structure */
    if (c->only_sites)
    {
        v = print_sites(b);
        if (v < 0)
        {
            return -1;
        }
    }
    else
    {
        v = pbwt_print(b, c->nohaps);
        if (v < 0)
        {
            return -1;
        }
    }

    /* Clean up allocated memory */
    pbwt_destroy(b);

    return 0;
}

int pbwt_print(const pbwt_t *b, const int nohaps)
{
    /* Check if pointer is NULL */
    if (b == NULL)
    {
        return -1;
    }

    size_t i = 0;
    size_t j = 0;

    for (i = 0; i < b->nsam; ++i)
    {
        /* Print sample identifier associated with haplotype i */
        if (b->sid[i])
        {
            if (nohaps)
            {
                printf("%.20s", b->sid[i]);
            }
            else
            {
                printf("%20.20s", b->sid[i]);
            }
        }
        else
        {
            fprintf(stderr, "pbwtmaster [ERROR]: problem reading sample identifier with index %5zu\n", i);
            return -1;
        }

        /* If a region is present */
        if (b->reg[i])
        {
            if (nohaps)
            {
                printf("\t%.30s", b->reg[i]);
            }
            else
            {
                printf("\t%30.30s", b->reg[i]);
            }
        }

        /* Print binary haplotype array for haplotype i */
        if (nohaps == 0)
        {
            putchar('\t');
            for (j = 0; j < b->nsite; ++j)
            {
                putchar((char)(b->data[TWODCORD(i, b->nsite, j)]));
            }
        }
        putchar('\n');
    }

    return 0;
}

int print_sites(const pbwt_t *b)
{
    size_t i = 0;

    for (i = 0; i < b->nsite; ++i)
    {
        printf("%s\t%s\t%lf\n", b->chr[i], b->rsid[i], b->cm[i]);
    }

    return 0;
}
