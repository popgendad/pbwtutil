#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ancmatch.h"

int pbwt_run(cmd_t *c)
{
    int v = 0;
    int no_query = 1;
    size_t i = 0;
    size_t nregs = 0;
    size_t nsites = 0;
    size_t qid = 0;
    khint_t k = 0;
    khash_t(floats) *result = NULL;
    char **reglist = NULL;
    pbwt_t *b = NULL;

    /* Initialize hash for results */
    result = kh_init(floats);
    v = 0;

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

    /* If no query name is given, pull index of query */
	for (i = 0; i < b->nsam; ++i)
	{
		if (strcmp(c->query, b->sid[i]) == 0)
		{
			qid = i;
            no_query = 0;
		}
        else
        {
            b->is_query[i] = 0;
        }
	}

    if (no_query)
    {
        fprintf(stderr, "ancmatch [ERROR]: cannot find haplotype with id %s\n", c->query);
        return 1;
    }

    /* Set query id */
    b->is_query[qid] = 1;

    /* Find all set-maximal matches */
    v = pbwt_query_match(b, c->minlen);
    match_regsearch(b, b->match, result, 0, nsites);

    /* Print hash */
    reglist = pbwt_get_reglist(b, &nregs);
    for (i = 0; i < nregs; ++i)
    {
        k = kh_get(floats, result, reglist[i]);
        if (kh_exist(result, k) && k != kh_end(result))
            fprintf(stdout, "%s\t%s\t%s\t%s\t%1.5lf\n", c->instub, b->sid[qid], b->reg[qid], reglist[i], kh_value(result, k));
        else
            fprintf(stdout, "%s\t%s\t%s\t%s\t0.00000\n", c->instub, b->sid[qid], b->reg[qid], reglist[i]);
    }

    free(reglist);

    /* Clean up allocated memory */
    pbwt_destroy(b);
    kh_destroy(floats, result);
    free(c->query);
    free(c);

    return v;
}
