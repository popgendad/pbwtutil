#include <stdio.h>
#include <stdlib.h>
#include "ancmatch.h"

int pbwt_print(const pbwt_t *, const int);

int pbwt_view(cmd_t *c)
{
    int v = 0;
    pbwt_t *b = NULL;

    b = pbwt_read(c->instub);
    v = pbwt_uncompress(b);
    v = pbwt_print(b, c->nohaps);
    pbwt_destroy(b);

    return v;
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
        size_t index = 0;

        index = b->ppa[i];

        /* Print prefix and divergence array member for haplotype i */
        printf("%5zu\t%5zu\t", b->div[i], index);

        /* Print binary haplotype array for haplotype i */
        if (nohaps == 0)
        {
	        for (j = 0; j < b->nsite; ++j)
	        {
	            putchar((char)(b->data[TWODCORD(index, b->nsite, j)]));
	        }
        }

        /* Print sample identifier associated with haplotype i */
        if (b->sid[index])
        {
            printf("\t%s", b->sid[index]);
        }
        else
        {
            fprintf(stderr, "ancmatch [ERROR]: problem reading sample identifier with index %5zu\n", index);
            return -1;
        }

        /* If a region is present */
        if (b->reg[index])
        {
            printf("\t%s\n", b->reg[index]);
        }
        else
        {
            putchar('\n');
        }
    }

    return 0;
}
