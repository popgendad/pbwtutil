#include <stdio.h>
#include <stdlib.h>
#include "pbwtmaster.h"

int pbwt_match(cmd_t *c)
{
    int v = 0;
    size_t i = 0;
    size_t nregs = 0;
    size_t qid = 0;
    khint_t k = 0;
    khash_t(floats) *result = NULL;
    khash_t(integer) *sdict = NULL;
    char **reglist = NULL;
    pbwt_t *b = NULL;

    /* Initialize hash for results */
    result = kh_init(floats);

    /* Read in the pbwt file from disk */
    b = pbwt_read(c->instub);
    if (b == NULL)
    {
        fprintf(stderr, "pbwtmaster [ERROR]: cannot read data from %s\n", c->instub);
        return -1;
    }

    /* Uncompress the haplotype data */
    pbwt_uncompress(b);

    /* Make dictionary of sample identifiers and their indices */
    sdict = pbwt_get_sampdict(b);
    if (sdict == NULL)
    {
        fputs("pbwtmaster [ERROR]: cannot construct sample identifier dictionary", stderr);
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

    /* Find all set-maximal matches */
    v = pbwt_query_match(b, c->minlen);
    if (v < 0)
    {
        fputs("pbwtmaster [ERROR]: error retrieving matches", stderr);
        return -1;
    }

    if (c->match_all)
    {
        match_print(b, b->match);
    }
    else
    {
        match_regsearch(b, b->match, result, 0, b->nsite);

        /* Get hash of regions */
        reglist = pbwt_get_reglist(b, &nregs);
        if (reglist == NULL)
        {
            fputs("pbwtmaster [ERROR]: cannot retrieve reg list", stderr);
            return -1;
        }

        /* Print region list to STDOUT */
        for (i = 0; i < nregs; ++i)
        {
            k = kh_get(floats, result, reglist[i]);
            if (kh_exist(result, k) && k != kh_end(result))
            {
                fprintf(stdout, "%s\t%s\t%s\t%s\t%1.5lf\n", c->instub, b->sid[qid], b->reg[qid], reglist[i], kh_value(result, k));
            }
            else
            {
                fprintf(stdout, "%s\t%s\t%s\t%s\t0.00000\n", c->instub, b->sid[qid], b->reg[qid], reglist[i]);
            }
        }
        free(reglist);
    }

    /* Clean up allocated memory */
    pbwt_destroy(b);
    kh_destroy(floats, result);
    kh_destroy(integer, sdict);
    free(c->query);
    free(c);

    return 0;
}
