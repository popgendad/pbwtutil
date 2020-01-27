#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pbwtmaster.h"

int pbwt_pileup(const cmd_t *c)
{
    int v = 0;
    size_t i = 0;
    size_t qid = 0;
    khint_t k = 0;
    khash_t(integer) *sdict = NULL;
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

    /* Make dictionary of sample identifiers and their indices */
    sdict = pbwt_get_sampdict(b);
    if (sdict == NULL)
    {
        fputs("pbwtmaster [ERROR]: cannot construct sample identifier dictionary\n", stderr);
        return -1;
    }

    /* Query user-input sample identifier */
    k = kh_get(integer, sdict, c->query);

    /* If sample ID is present, mark haplotype as query */
    if (kh_exist(sdict, k) && k != kh_end(sdict))
    {
        qid = kh_value(sdict, k);
        b->is_query[qid] = TRUE;
    }
    else
    {
        fprintf(stderr, "pbwtmaster [ERROR]: cannot find haplotype with id %s\n", c->query);
        return -1;
    }

    v = pbwt_all_query_match(b, c->minlen, add_interval);
    if (b->intree == NULL)
    {
    	puts("Problem with reporting to interval tree");
    	return -1;
    }

    for (i = 0; i < b->nsite - 10; i += 10)
    {
        size_t start_pos = i;
        size_t end_pos = i + 10;
        size_t c = 0;
        match_count(b, b->intree, &c, start_pos, end_pos);
        printf("%zu\t%zu\t%zu\n", start_pos, end_pos, c);
    }

	kh_destroy(integer, sdict);
    pbwt_destroy(b);

	return 0;
}
