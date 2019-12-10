#include <stdio.h>
#include <stdlib.h>
#include "pbwtmaster.h"

int pbwt_print(const pbwt_t *, const int);
int print_sites(const pbwt_t *);

int pbwt_view(cmd_t *c)
{
    int v = 0;
    pbwt_t *b = NULL;

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
        fputs("pbwtmaster [ERROR]: error uncompressing haplotype data", stderr);
        return -1;
    }

    v = pbwt_build(b);

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
    free(c);

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
        size_t index = b->ppa[i];

        /* Print sample identifier associated with haplotype i */
        if (b->sid[index])
        {
            printf("%s", b->sid[index]);
        }
        else
        {
            fprintf(stderr, "pbwtmaster [ERROR]: problem reading sample identifier with index %5zu\n", index);
            return -1;
        }

        /* If a region is present */
        if (b->reg[index])
        {
            printf("\t%s", b->reg[index]);
        }
        putchar('\t');

        /* Print binary haplotype array for haplotype i */
        if (nohaps == 0)
        {
            for (j = 0; j < b->nsite; ++j)
            {
                putchar((char)(b->data[TWODCORD(index, b->nsite, j)]));
            }
        }

        /* Print prefix and divergence array member for haplotype i */
        printf("%5zu\t%5zu\n", b->div[i], index);
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
